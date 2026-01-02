function [cur_rts,cur_intens,cur_mono_isointens] = get_histone2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,hno,special)
%%
% GET_HISTONE2  Two-component deconvolution (2 isoforms) of a histone peptide using MS1+MS2 evidence.
%
% Synopsis:
%   [cur_rts, cur_intens, cur_mono_isointens] = get_histone2( ...
%       MS1_index, MS1_peaks, MS2_index, MS2_peaks, ptol, unitdiff, Mods, His, hno, special)
%
% High-level purpose:
%   This routine resolves a target peptide signal into TWO co-eluting/partially overlapping components
%   (e.g., two modification isomers "iso-1" and "iso-2"). It:
%     1) Quantifies the *total* MS1 area for the peptide (all charge states) via GET_HISTONE1.
%     2) Estimates the *relative ratio* between the two components using MS2 key fragments
%        (MatchMS2 → key ion sets K1/K2) within a narrow RT window and/or a smooth fallback.
%     3) Splits the total intensities across the two components using that ratio.
%   Output intensities are reported as a 2×nCharge matrix (row 1 for iso-1, row 2 for iso-2).
%
% Key behavior:
%   - If the base peptide (from GET_HISTONE1) has zero intensity at charge 1, the function returns early.
%   - For the 2-component ratio, the helper GET_RATIO_2ISO uses MS2 evidence when available; otherwise
%     it builds a smooth (cosine) heuristic profile across the RT segment to keep numerical stability.
%   - cur_mono_isointens returns TWO mono-XIC-like traces (num_MS1×2) corresponding to the split signals.
%
% Inputs (brief):
%   MS1_index, MS1_peaks : MS1 scan index and peak table (see codebase conventions).
%   MS2_index, MS2_peaks : MS2 scan index and peak table (with precursor m/z, RT, instrument code, etc.).
%   ptol                 : ppm tolerance (100 is coerced to 10 as per historical convention).
%   unitdiff             : 13C spacing in Da (~1.003355), used as unitdiff / charge for isotopic seeds.
%   Mods                 : modification catalog (passed to get_key_ions*/K1,K2 selection).
%   His                  : struct with peptide info (pep_mz, pep_ch, rt_ref, mod_short, outpath, outfile...).
%   hno                  : row index of the FIRST isoform (iso-1); code expects paired isoforms (hno,hno+1).
%   special              : control flags (e.g., nDAmode: 1=DDA, 2=DIA; nsource/nsubtype affect NH handling).
%
% Outputs:
%   cur_rts             : [2 × nCharge] apex RTs for iso-1 and iso-2 (rows 1 and 2), repeated across charges.
%   cur_intens          : [2 × nCharge] areas for iso-1 and iso-2 per charge (split from total via ratio).
%   cur_mono_isointens  : [num_MS1 × 2] split mono-isotopic profiles (column 1 for iso-1, column 2 for iso-2).
%
% Calls:
%   get_histone1 → baseline MS1 quant; get_ratio_2iso → MS2-informed two-component split; MatchMS2 (nested).
%

[npep,ncharge] = size(His.pep_mz);%#ok
cur_rts = zeros([2,ncharge]);
cur_intens = zeros([2,ncharge]);
num_MS1 = size(MS1_index,1);
cur_mono_isointens = zeros([num_MS1,2]);

% 1) Get total MS1 signal (per charge) around the reference RT using get_histone1.
[h_rts,h_intens] = get_histone1(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
if 0==h_intens(1)
    % If the primary charge state has zero area, bail out (nothing to split).
    return;
end;

% 2) Estimate the two-component ratio (and their RTs) using MS2 evidence (or a smooth fallback).
[s_rts,ratio,cur_mono_isointens] = get_ratio_2iso(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,hno,h_rts,special);

% 3) Broadcast the two RTs across all charge states, and split total intensity by the ratio.
cur_rts(1,1:ncharge) = repmat(s_rts(1),[1,ncharge]);      % iso-1 apex RT (replicated across charges)
cur_rts(2,1:ncharge) = repmat(s_rts(2),[1,ncharge]);      % iso-2 apex RT
cur_intens(1,1:ncharge) = h_intens*ratio;                 % iso-1 shares 'ratio' of the total
cur_intens(2,1:ncharge) = h_intens*(1-ratio);             % iso-2 gets the complement

function [s_rts,ratio,cur_mono_isointens] = get_ratio_2iso(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,hno,h_rts,special)
%%
% GET_RATIO_2ISO  Build a two-component time profile and ratio using MS2 key-ion evidence (or a heuristic).
%
% Returns:
%   s_rts  : [1×2] apex RTs for iso-1 and iso-2.
%   ratio  : scalar in [0,1], area fraction assigned to iso-1 (iso-2 gets 1−ratio).
%   cur_mono_isointens : [num_MS1×2] split mono-XICs, aligned to MS1 scans and filled with zeros outside IX.
%
% Strategy:
%   • Extract an MS1 mono-XIC for the first charge state within ±0.5 min around rt_ref(hno).
%   • Determine the RT segment [rt1,rt2] (DDA: full segment; DIA: aligned to the target isolation window).
%   • Obtain per-MS2-scan ratios from key fragment sets (K1 for iso-1, K2 for iso-2) via MatchMS2.
%   • Map those ratios to the MS1 scan grid; if sparse, use a smooth cosine fallback with small amplitude.
%   • Split the mono-XIC into two curves (n_inten2, n_inten3) with ratio1 and (1−ratio1), smooth, zero low tails.
%   • Pick apex times of both curves and integrate areas to compute the global split 'ratio'.

npart = 2;
s_rts = zeros([1,npart]);
ratio = 0;
num_MS1 = size(MS1_index,1);
cur_mono_isointens = zeros([num_MS1,npart]);

% MS1 extraction window (tight) around rt_ref(hno): ±0.5 min.
delta = 0.5;
p = find( MS1_index(:,2)>=His.rt_ref(hno)-delta );
rt_i1 = p(1);
pp = find( MS1_index(:,2)<=His.rt_ref(hno)+delta );
rt_i2 = pp(end);

% Build mono-XIC for the first charge state using 4 isotope seeds {M−1, M, M+1, M+2}.
c_mz = His.pep_mz(hno,1);
c_ch = His.pep_ch(hno,1);
c_ref_isomzs = [c_mz-unitdiff/c_ch c_mz c_mz+unitdiff/c_ch c_mz+2*unitdiff/c_ch];
if ptol==100
    ptol = 10;
end;
if ptol>100 && c_ch>=3
    nC13 = 1;
else
    nC13 = 0;
end;
[c_isorts,c_ref_isointens] = GetProfiles(MS1_index,MS1_peaks,c_ref_isomzs,c_ch,ptol,nC13,rt_i1:rt_i2);
j = 2;
c_mono_isointens = c_ref_isointens(:,j);

% Define the RT segment [rt1, rt2] inside this window (later refined for DIA).
IX = rt_i1:rt_i2;
% (Original code had a local-maximum finder commented out; we keep the simple end-to-end segment.)
rt1 = c_isorts(IX(1));
rt2 = c_isorts(IX(end));

% DIA specialization: restrict [rt1, rt2] to the subsequence actually covered by the target MS2 isolation.
if 2==special.nDAmode
    premzs = unique(MS2_index(:,4));
    [tmp,ii] = min( abs(premzs-c_mz) );%#ok
    target = premzs(ii);
    num_MS2 = size(MS2_index,1);
    flag = zeros([num_MS2,1]);
    p = find( MS2_index(:,2)>=rt1 );
    if 1==isempty(p)
        return;
    end;
    i1 = p(1);
    pp = find( MS2_index(:,2)<=rt2 );
    if 1==isempty(pp)
        return;
    end;
    i2 = pp(end);
    for i=i1:i2
        cen_mz = MS2_index(i,4);
        if 0==cen_mz-target
            flag(i) = 1;
        end;
    end;
    
    % Map DIA MS2 scans to MS1 scans to refine the MS1 time bounds.
    ms2pos = find(flag==1);
    if 1==isempty(ms2pos)
        return;
    end;
    IX = zeros([length(ms2pos),1]);
    for i=1:length(ms2pos)
        IX(i) = find(MS1_index(:,1)==MS2_index(ms2pos(i),1));
    end;
    rt1 = c_isorts(IX(1));
    rt2 = c_isorts(IX(end));
end;

% Compute per-MS2-scan ratios using diagnostic fragment sets for iso-1 vs iso-2.
[ms2pos,ms2ratios] = MatchMS2(MS2_index,MS2_peaks,c_mz,c_ch,ptol,unitdiff,Mods,His,hno,rt1,rt2,special);

% Build a time-resolved ratio vector 'ratio1' aligned to MS1 scans:
nlen = length(IX);
ratio1 = zeros([nlen,1]);

if length(ms2pos)<nlen/2
    % Sparse MS2 evidence → use a smooth, small-amplitude cosine profile as a stabilizing fallback.
    if nlen<=3
        return;
    end;
    for i=1:nlen
        ratio1(i) = cos((i-1)*pi/(2*nlen)); % from 1 → 0 smoothly across the window
    end;
    % Directionality: decide which iso peaks first using rt_ref ordering of (hno) vs (hno+1).
    if His.rt_ref(hno)>His.rt_ref(hno+1)
        ratio1 = rot90(rot90(ratio1)); % reverse vector
    end;
    ratio1 = ratio1/10;                % keep amplitude low (heuristic)
    ratio2 = 1-ratio1;
else
    % Map MS2 ratios to the MS1 scan grid (by scan number intersection).
    ms1scans = MS2_index(ms2pos,1);
    for i=1:nlen
        c_scan = MS1_index(IX(i),1);
        [tf,loc] = ismember(c_scan,ms1scans);
        if 1==tf
            ratio1(i) = ms2ratios(loc);
        end;
    end;

    % Trim low-ratio tails beyond the main maximum to encourage bimodality/separation.
    if His.rt_ref(hno)<=His.rt_ref(hno+1)
        [tmp,x1] = max(ratio1);
        nX1 = x1:nlen;
        x2 = find(ratio1(nX1)<0.1);
        if tmp>0.1 && 0==isempty(x2)
            nX2 = nX1(x2(1)):nlen;
            ratio1(nX2) = 0.001;
        end;
    else
        [tmp,x1] = max(ratio1);
        nX1 = 1:x1;
        x2 = find(ratio1(nX1)<0.1);
        if tmp>0.1 && 0==isempty(x2)
            nX2 = 1:nX1(x2(end));
            ratio1(nX2) = 0.001;
        end;
    end;
    ratio1 = smooth(ratio1,3);
    ratio2 = 1-ratio1;
end;

% Split the mono-XIC by the ratio and smooth the two resulting curves.
n_rt = c_isorts(IX);
n_inten1 = c_mono_isointens(IX);
n_inten2 = n_inten1.*ratio1;           % iso-1 time-resolved intensity
n_inten3 = n_inten1.*ratio2;           % iso-2 time-resolved intensity

new_inten2 = smooth(n_inten2,3);
new_inten3 = smooth(n_inten3,3);
new_inten2(find(new_inten2<0.2*max(new_inten2))) = 0;%#ok
new_inten3(find(new_inten3<0.2*max(new_inten3))) = 0;%#ok

% Pick apex positions for both curves (using GetTopBottom11 on the smoothed signals).
[nt2,nb2,top1_idx2] = GetTopBottom11(new_inten2);%#ok
[nt3,nb3,top1_idx3] = GetTopBottom11(new_inten3);%#ok
if 1==isempty(nt2)
    top1 = His.rt_ref(hno);
else
    top1 = n_rt(nt2(top1_idx2));
end;
if 1==isempty(nt3)
    top2 = His.rt_ref(hno);
else
    top2 = n_rt(nt3(top1_idx3));
end;

% Integrate areas of the two components over the segment via spline sampling at 0.005 min.
xx = n_rt(1):0.005:n_rt(end);
yy = spline(n_rt,n_inten2,xx);
area1 = sum(abs(yy));

xx = n_rt(1):0.005:n_rt(end);
yy = spline(n_rt,n_inten3,xx);
area2 = sum(abs(yy));

s_rts = [top1 top2];
ratio = area1/(eps+area1+area2);     % global split assigned to iso-1 (iso-2 gets 1−ratio)

% Build full-length (num_MS1) split mono-XICs, zero-padded outside [IX(1):IX(end)].
if 2~=special.nDAmode
    cur_mono_isointens(:,1) = [zeros([IX(1)-1,1]);n_inten2;zeros([num_MS1-IX(end),1])];
    cur_mono_isointens(:,2) = [zeros([IX(1)-1,1]);n_inten3;zeros([num_MS1-IX(end),1])];
else
    % In DIA, resample the split curves back onto the MS1 time grid for the IX subrange.
    xx = c_isorts(IX(1):IX(end));
    yy = spline(n_rt,n_inten2,xx);
    cur_mono_isointens(:,1) = [zeros([IX(1)-1,1]);yy;zeros([num_MS1-IX(end),1])];
    yy = spline(n_rt,n_inten3,xx);
    cur_mono_isointens(:,2) = [zeros([IX(1)-1,1]);yy;zeros([num_MS1-IX(end),1])];
end;

% Produce a small PDF diagnostic plot with the experimental mono-XIC and its two split components.
set(gcf,'visible','off');
out_file1 = fullfile(His.outpath,['Iso_',His.outfile,'_',His.mod_short{hno},'_',His.mod_short{hno+1},'.pdf']);
plot(n_rt,n_inten1,'linestyle','-','linewidth',2,'color','k');
hold on;
plot(n_rt,n_inten2,'linestyle','--','linewidth',2,'color','r');
plot(n_rt,n_inten3,'linestyle','-.','linewidth',2,'color','b');
xlabel('time (min)');
ylabel('intensity');
legend('experiment',His.mod_short{hno},His.mod_short{hno+1});
if area1/(eps+area1+area2)<0.03
    r1 = floor(1000*area1/(eps+area1+area2))*0.1;
else
    r1 = floor(100*area1/(eps+area1+area2));
end;
r2 = 100-r1;
title([His.mod_short{hno},'/',His.mod_short{hno+1},':',num2str(r1),'%:',num2str(r2),'%']);
print('-dpdf',out_file1);
close();
%}

function [ms2pos,ms2ratios] = MatchMS2(MS2_index,MS2_peaks,c_mz,c_ch,ptol,unitdiff,Mods,His,hno,rt1,rt2,special)
%%
% MATCHMS2  Identify relevant MS2 scans and compute per-scan iso-1 vs iso-2 ratios from key ions.
%
% Returns:
%   ms2pos    : indices of MS2 scans considered relevant within [rt1,rt2].
%   ms2ratios : per-scan ratios in [0,1], computed from key ion sets (K1 for iso-1, K2 for iso-2).
%
% Modes:
%   - DDA (special.nDAmode==1): accept MS2 scans whose precursor m/z matches any of {M, M+1, M+2}/z.
%   - DIA (special.nDAmode==2): accept MS2 scans whose isolation center equals the nearest DIA window to c_mz.
%   - Else: return empty (no MS2 support).
%
% Instrument handling:
%   MS2_index(:,6) encodes instrument types (IT/FT; CID/ETD/HCD).
%   Odd codes (…IT) → tol = 0.4 Da; even codes (…FT) → tol = 0.02 Da.
%   ActiveType = 'ETD' if code in {3,4}; otherwise 'CID' (covers CID/HCD for key-ion logic here).
%
% Key ions:
%   - get_key_ions / get_key_ionsH (NH-mass variant) produce matched lists K1 and K2 (same length),
%     diagnostic for iso-1 vs iso-2 given Mods, sequences, and ActiveType.

% NH-mass flag (some sources/subtypes use a different protonation mass handling).
if 4==special.nsource && (0~=special.nsubtype && 2~=special.nsubtype)
    nhmass = 1;
else
    nhmass = 0;
end;

% --- Select MS2 scans by mode ---
num_MS2 = size(MS2_index,1);
if 1==special.nDAmode
    % DDA: allow precursor at M, M+1, M+2 (accounting for 13C isotopes) within ppm tolerance.
    sets = [0 1 2];
    mzs = c_mz + sets*unitdiff/c_ch;
    flag = zeros([num_MS2,1]);
    p = find( MS2_index(:,2)>=rt1 );
    if 1==isempty(p)
        ms2pos = [];
        ms2ratios = [];
        return;
    end;
    i1 = p(1);
    pp = find( MS2_index(:,2)<=rt2 );
    if 1==isempty(pp)
        ms2pos = [];
        ms2ratios = [];
        return;
    end;
    i2 = pp(end);
    for i=i1:i2
        cen_mz = MS2_index(i,4);
        ix = find(abs(mzs-cen_mz)<ptol*cen_mz*1e-6);%#ok
        if 0==isempty(ix)
            flag(i) = 1;
        end;
    end;

    % Gather selected MS2 scans.
    ms2pos = find(flag==1);
    if 1==isempty(ms2pos)
        ms2ratios = [];
        return;
    end;
elseif 2==special.nDAmode
    % DIA: match the DIA isolation window whose center equals the nearest to c_mz.
    premzs = unique(MS2_index(:,4));
    [tmp,ii] = min( abs(premzs-c_mz) );%#ok
    target = premzs(ii);
    flag = zeros([num_MS2,1]);
    p = find( MS2_index(:,2)>=rt1 );
    if 1==isempty(p)
        ms2pos = [];
        ms2ratios = [];
        return;
    end;
    i1 = p(1);
    pp = find( MS2_index(:,2)<=rt2 );
    if 1==isempty(pp)
        ms2pos = [];
        ms2ratios = [];
        return;
    end;
    i2 = pp(end);
    for i=i1:i2
        cen_mz = MS2_index(i,4);
        if 0==cen_mz-target
            flag(i) = 1;
        end;
    end;

    % Gather selected MS2 scans.
    ms2pos = find(flag==1);
    if 1==isempty(ms2pos)
        ms2ratios = [];
        return;
    end;
else
    ms2pos = [];
    ms2ratios = [];
    return;
end;

% --- Determine instrument behavior and key ions K1/K2 ---
instruments = MS2_index(ms2pos,6);% MS2dirs = {'CIDIT','CIDFT','ETDIT','ETDFT','HCDIT','HCDFT'};
if 1==length(unique(instruments))
    % Uniform instrument across scans: define ActiveType and m/z tolerance once.
    c_instrument = instruments(1);
    if 3==c_instrument || 4==c_instrument
        ActiveType = 'ETD';
    else
        ActiveType = 'CID';
    end;
    if 1==mod(c_instrument,2)
        tol = 0.4;
    else
        tol = 0.02;
    end;

    % Key ions for iso-discrimination:
    if 1==nhmass
        [K1,K2] = get_key_ionsH(His,hno,hno+1,Mods,ActiveType);
    else
        [K1,K2] = get_key_ions(His,hno,hno+1,Mods,ActiveType);
    end;
end;

% Precompute spectrum row index starts for each MS2 scan.
index = [1;MS2_index(1:num_MS2,7)];
ms2ratios = zeros([1,length(ms2pos)]);
for i=1:length(ms2pos)
    cno = ms2pos(i);
%   In DIA you could also consider neighboring scans (cno-1:cno+1) – original alt is commented out.
    newpos = cno;
    for pno = newpos
        if pno<1 || pno>num_MS2
            continue;
        end;
        if 1<length(unique(instruments))
            % Mixed instruments: update ActiveType/tol per scan.
            c_instrument = MS2_index(pno,6);
            if 3==c_instrument || 4==c_instrument
                ActiveType = 'ETD';
            else
                ActiveType = 'CID';
            end;
            if 1==mod(c_instrument,2)
                tol = 0.4;
            else
                tol = 0.02;
            end;
        end;

        if 1<length(unique(instruments))
            % If instruments mixed, recompute key ions accordingly.
            if 1==nhmass
                [K1,K2] = get_key_ionsH(His,hno,hno+1,Mods,ActiveType);
            else
                [K1,K2] = get_key_ions(His,hno,hno+1,Mods,ActiveType);
            end;
        end;

        % Read spectrum and find intensities at each key ion with ±tol matching.
        IX = index(pno):index(pno+1)-1;
        mz = MS2_peaks(IX,1);
        inten = MS2_peaks(IX,2);

        intens1 = zeros([1,length(K1)]);
        intens2 = zeros([1,length(K1)]);
        for j=1:length(K1)
            ix1 = find(abs(mz-K1(j))<=tol);
            ix2 = find(abs(mz-K2(j))<=tol);
            [tmp,x1] = min(abs(mz(ix1)-K1(j)));%#ok
            [tmp,x2] = min(abs(mz(ix2)-K2(j)));%#ok
            % (Alternative max-intensity matching is commented; here we pick the nearest m/z.)
            if 0==isempty(ix1)
                intens1(j) = inten(ix1(x1));
            end;
            if 0==isempty(ix2)
                intens2(j) = inten(ix2(x2));
            end;
        end;
        % Per-scan ratio: iso-1 contribution over (iso-1 + iso-2), robustified by sum over key ions.
        ratio1 = sum(intens1)/(eps+sum(intens1)+sum(intens2));
        if ms2ratios(i)<ratio1
            ms2ratios(i) = ratio1;
        end;
    end;
end;
