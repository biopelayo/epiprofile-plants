function [rts,top1_rt,inten_sum,top1_inten_sum] = get_rts22( ...
    MS1_index, MS1_peaks, MS2_index, MS2_peaks, ...
    ptol, unitdiff, His, hno, nsplit, t1, t2, nhmass)
%% get_rts22 — MS1+MS2 co-evidence peak picking with “near-miss” rescue of top1
%
% GOAL (plain English, non-expert friendly)
%   Given a target histone peptide (indexed by 'hno' in 'His'), this function
%   finds its likely elution time(s) within a retention-time window [t1, t2].
%   It:
%     1) extracts the peptide’s MS1 monoisotopic XIC (using isotopic channels),
%     2) segments the window into candidate elution features and quantifies them,
%     3) checks for diagnostic MS2 fragments around each candidate,
%     4) scores MS1–MS2 shape similarity and selects a best (top-1) apex.
%
%   Compared to get_rts2, get_rts22 adds a “near-miss rescue”:
%     if the originally selected apex is very close to the max MS2 evidence
%     (within 1–2 similarity units), it can be promoted to the winning bin.
%
% INPUTS (minimal expectations; same conventions as in EpiProfile2)
%   MS1_index: [N_MS1 x ≥2]
%       (1) MS1 scan id, (2) RT (same units as t1,t2). Other cols ignored.
%   MS1_peaks: [Σ_MS1_peaks x 2]  (m/z, intensity) concatenated over scans.
%   MS2_index: [N_MS2 x ≥8]
%       (1) MS2 scan id, (2) RT, (4) precursor m/z, (6) instrument code,
%       (7) row pointer into MS2_peaks, (8) noise/threshold estimate.
%   MS2_peaks: [Σ_MS2_peaks x 2]  (m/z, intensity) concatenated over scans.
%   ptol:      MS1 tolerance (ppm or special code; ptol==100 → internally 10).
%   unitdiff:  ~1.00335 Da (13C mass difference) to position isotope channels.
%   His:       struct with His.pep_mz(:,1) and His.pep_ch(:,1) at least.
%   hno:       peptide index in 'His'.
%   nsplit:    segmentation mode (both branches call same helper here).
%   t1,t2:     RT window to search (must overlap the run RTs).
%   nhmass:    0/1 toggling which key-ion generator to use.
%
% OUTPUTS
%   rts            : candidate apex RTs (one per detected feature).
%   top1_rt        : best RT after combining MS1 and MS2 evidence.
%   inten_sum      : MS1 area-like intensity per candidate (same length as rts).
%   top1_inten_sum : MS1 sum for the chosen top-1 (penalized /1e4 if weak MS2).
%
% DEPENDENCIES (not defined here; must exist on MATLAB path)
%   GetProfiles, GetTopBottom11, GetMods, get_key_ions1, get_key_ions1H
%
% CAVEATS
%   - MS2 precursor matching below is by exact equality to a chosen target m/z.
%     If your data has rounding variation, change that to a ppm/Da tolerance.
%   - RT units: ensure MS1, MS2, t1, t2 all share the same time unit.
%   - If top1_rt < 4 (very early), returns empty (solvent-front heuristic).

%% 1) Window validity check (fast exit if invalid)
num_MS1 = size(MS1_index,1);
if t1>MS1_index(num_MS1,2) || t2<0 || t2<t1
    rts = []; top1_rt = []; inten_sum = []; top1_inten_sum = [];
    return;
end

%% 2) Resolve target peptide and define isotope channels
% Monoisotopic m/z and charge from the histone library
c_mz = His.pep_mz(hno,1);
c_ch = His.pep_ch(hno,1);

% Convert RT limits into MS1 row indices for the slice
p     = find(MS1_index(:,2)>=t1);  rt_i1 = p(1);
pp    = find(MS1_index(:,2)<=t2);  rt_i2 = pp(end);

% Reference isotope m/zs: M−1, M, M+1, M+2 (spaced by unitdiff/charge)
c_ref_isomzs = [c_mz-unitdiff/c_ch, c_mz, c_mz+unitdiff/c_ch, c_mz+2*unitdiff/c_ch];

% Normalize ptol conventions (legacy EpiProfile behavior)
if ptol==100, ptol = 10; end
if ptol>100 && c_ch>=3
    nC13 = 1;  % enable 13C handling for highly charged precursors under special ptol
else
    nC13 = 0;
end

%% 3) Extract MS1 profiles and build monoisotopic trace
% GetProfiles returns: c_isorts (RT vector), c_ref_isointens (nScans x 4)
[c_isorts,c_ref_isointens] = GetProfiles( ...
    MS1_index, MS1_peaks, c_ref_isomzs, c_ch, ptol, nC13, rt_i1:rt_i2);

% Use monoisotopic channel (column 2) as main XIC
c_mono_isointens = c_ref_isointens(:,2);

%% 4) Segment the trace into peaks and compute intensities
% GetTopBottom11 returns: nt (apex indices), nb (boundaries), top1_idx, inten_sum
if 1==nsplit
    [nt,nb,top1_idx,inten_sum] = GetTopBottom11(c_mono_isointens);
else
    [nt,nb,top1_idx,inten_sum] = GetTopBottom11(c_mono_isointens);
end
rts = c_isorts(nt);  % candidate RTs

%% 5) Combine with MS2 evidence (diagnostic fragments + cosine similarity)
Mods = GetMods();        % modification catalog (used to compute key ions)
flen = 0;                % tracks number of expected neutral-type ion families

% Optional denoising: zero out very small intensities (<1.5% max)
% (Original code kept smoothing commented; we preserve that choice.)
% w = smooth(c_mono_isointens,3);
w = c_mono_isointens;
w(find(w<0.015*max(w))) = 0; %#ok<FNDSB>

similarity = zeros([length(nt),1]);  % per-candidate similarity accumulator

for ino = 1:length(nt)
    % Boundaries for this candidate region
    i1 = nb(ino); i2 = nb(ino+1);

    % Trim to first/last non-zero within the region to avoid flat tails
    while i1<=length(w) && 0==w(i1), i1 = i1+1; end
    if i1>nb(ino),           i1 = i1-1; end
    while i2>=1 && 0==w(i2), i2 = i2-1; end
    if i2<nb(ino+1),         i2 = i2+1; end

    IX  = i1:i2;
    rt1 = c_isorts(IX(1));     % local region start RT
    rt2 = c_isorts(IX(end));   % local region end   RT

    % Pull MS2 scans & intensities for diagnostic ions matching this peptide
    [ms2pos,ms2rts,ms2intens,posn,posc] = ...
        MatchMS2(MS2_index,MS2_peaks,Mods,His,hno,rt1,rt2,nhmass);

    % Track how many neutral-type ion families exist (used later for gating)
    if flen < length(posn), flen = length(posn); end

    % If no MS2 support around this region, skip similarity scoring
    if isempty(ms2pos), continue; end

    % Define a symmetric RT window around the apex based on median MS2 spacing
    nstep = median(ms2rts(ms2pos(2:end)) - ms2rts(ms2pos(1:end-1))) * 3.5;
    nst   = c_isorts(nt(ino)) - nstep;
    ntm   = c_isorts(nt(ino)) + nstep;

    % Subselect MS2 scans within that window
    np1 = find(ms2rts(ms2pos) >= nst);
    np2 = find(ms2rts(ms2pos) <= ntm);
    X   = np1(1):np2(end);

    % Map those MS2 scans to the nearest MS1 scan indices to sample the XIC
    XX = zeros([1,length(X)]);
    for jno=1:length(X)
        XX(jno) = find(MS1_index(:,1) == MS2_index(ms2pos(X(jno)),1));
    end

    % Pull the MS1 monoisotopic intensities at aligned indices and smooth lightly
    p_intens = c_mono_isointens(XX);
    p_intens = smooth(p_intens,3);

    % Build a conservative similarity threshold and add tiny baseline to avoid zeros
    [sim_thr,ex_intens] = get_thr(p_intens);
    p_intens = p_intens + ex_intens;

    % Compare MS1 shape vs each neutral-type diagnostic ion family
    for kno=1:length(posn)
        f_intens = ms2intens(ms2pos(X),kno);
        f_intens = smooth(f_intens,3) + ex_intens;
        e_sim = sum(p_intens.*f_intens) / sqrt( sum(p_intens.*p_intens) * sum(f_intens.*f_intens) );
        if e_sim > sim_thr
            similarity(ino) = similarity(ino) + e_sim;
        end
    end

    % Do the same for complementary/charge-dependent ions (posc)
    for kno=1:length(posc)
        qno     = kno + length(posn);
        f_intens = ms2intens(ms2pos(X),qno);
        f_intens = smooth(f_intens,3) + ex_intens;
        e_sim = sum(p_intens.*f_intens) / sqrt( sum(p_intens.*p_intens) * sum(f_intens.*f_intens) );
        if e_sim > sim_thr
            similarity(ino) = similarity(ino) + e_sim;
        end
    end
end

%% 6) Select the winning apex (top-1) — with “near-miss rescue”
% ‘tmp’ is the integer ceiling of the best similarity observed across candidates.
tmp = max(ceil(similarity));

if tmp > flen/2
    % NEW vs get_rts2: if the apex chosen by MS1 segmentation (top1_idx)
    % is within 1–2 similarity “units” of the maximum, promote it to the max.
    if (ceil(similarity(top1_idx))==tmp-1 && tmp>1) || ...
       (ceil(similarity(top1_idx))==tmp-2 && tmp>2)
        similarity(top1_idx) = tmp;
    end

    % Among all candidates that reach that max similarity, pick the one
    % with the largest MS1 integrated intensity (inten_sum).
    ii = find(ceil(similarity)==tmp);
    [~,id] = max(inten_sum(ii)); %#ok<ASGLU>
    top1_idx = ii(id);

    top1_rt        = c_isorts(nt(top1_idx));
    top1_inten_sum = inten_sum(top1_idx);
else
    % If overall MS2 support is weak, keep the original top1 but penalize intensity
    top1_rt        = c_isorts(nt(top1_idx));
    top1_inten_sum = inten_sum(top1_idx) / 1e4;
end

%% 7) Solvent-front guard (very early RTs are unreliable)
if top1_rt < 4
    rts = []; top1_rt = []; inten_sum = []; top1_inten_sum = [];
    return;
end

end % ===== end of main get_rts22 =====


%% ------------------ Local helper: dynamic similarity threshold ------------------
function [sim_thr,ex_intens] = get_thr(p_intens0)
% PURPOSE
%   Build a conservative (local) similarity threshold by comparing against three
%   tiny baseline patterns (flat / rising / falling). The minimum of the three
%   cosine values is used as a threshold; we also return the baseline vector
%   added to both MS1/MS2 shapes to avoid zero-division and spurious cosines.

% Three epsilon-level patterns
x(1:length(p_intens0),1) = (ones(1,length(p_intens0)))'*eps;
x(1:length(p_intens0),2) = (1:length(p_intens0))'*eps;
x(1:length(p_intens0),3) = (length(p_intens0):-1:1)'*eps;

% Cosine vs flat
p_intens = p_intens0 + x(:,1);
f_intens = x(:,1);
sim_thrs(1) = sum(p_intens.*f_intens)/sqrt( sum(p_intens.*p_intens)*sum(f_intens.*f_intens) );

% Cosine vs rising
p_intens = p_intens0 + x(:,2);
f_intens = x(:,2);
sim_thrs(2) = sum(p_intens.*f_intens)/sqrt( sum(p_intens.*p_intens)*sum(f_intens.*f_intens) );

% Cosine vs falling
p_intens = p_intens0 + x(:,3);
f_intens = x(:,3);
sim_thrs(3) = sum(p_intens.*f_intens)/sqrt( sum(p_intens.*p_intens)*sum(f_intens.*f_intens) );

% Final threshold and baseline used
[sim_thr,i] = min(sim_thrs);
ex_intens   = x(:,i);
end


%% ------------------ Local helper: gather MS2 key-ion intensities ----------------
function [ms2pos,ms2rts,ms2intens,posn,posc,ActiveType,K1] = ...
    MatchMS2(MS2_index,MS2_peaks,Mods,His,hno,rt1,rt2,nhmass)
% WHAT IT DOES
%   1) Selects MS2 scans whose precursor m/z equals the target’s (closest unique).
%   2) Keeps those within [rt1, rt2].
%   3) Infers fragmentation type from instrument code, sets m/z tolerance.
%   4) Builds diagnostic ions (K1) using get_key_ions1/1H and extracts their intensities.

% --- Select scans in the local RT window
num_MS2 = size(MS2_index,1);
c_mz = His.pep_mz(hno,1);

% Choose the closest precursor m/z present in the run as the “target”
premzs = unique(MS2_index(:,4));
[~,ii] = min( abs(premzs-c_mz) ); %#ok<ASGLU>
target = premzs(ii);

% Guard for empty local window
flag = zeros([num_MS2,1]);
p  = find( MS2_index(:,2)>=rt1 );
pp = find( MS2_index(:,2)<=rt2 );
if isempty(p) || isempty(pp)
    ms2pos=[]; ms2rts=[]; ms2intens=[]; posn=[]; posc=[]; ActiveType=[]; K1=[];
    return;
end
i1 = p(1); i2 = pp(end);

% Exact-m/z match to target inside the window (change to ppm/Da if needed)
for i=i1:i2
    cen_mz = MS2_index(i,4);
    if 0==cen_mz-target
        flag(i) = 1;
    end
end
ms2pos = find(flag==1);

% Require >4 scans to proceed (empirical stability)
if length(ms2pos)<=4
    ms2pos=[]; ms2rts=[]; ms2intens=[]; posn=[]; posc=[]; ActiveType=[]; K1=[];
    return;
end

% --- Fragmentation type and tolerance from instrument code
ms2rts = MS2_index(:,2);
instruments = MS2_index(ms2pos,6); % odd≈trap, even≈FT; 3/4≈ETD family

c_instrument = instruments(1);
if 3==c_instrument || 4==c_instrument
    ActiveType = 'ETD';
else
    ActiveType = 'CID';
end
if mod(c_instrument,2)==1
    tol = 0.4;   % Da, trap-like
else
    tol = 0.02;  % Da, FT-like
end

% --- Build diagnostic ion list
if 1==nhmass
    [K1,posn,posc] = get_key_ions1H(His,hno,Mods,ActiveType);
else
    [K1,posn,posc] = get_key_ions1(His,hno,Mods,ActiveType);
end

% --- Extract intensities for those ions per selected scan
index = [1; MS2_index(1:num_MS2,7)];           % row pointer into MS2_peaks
ms2intens = zeros([num_MS2,length(K1)]);

for i=1:length(ms2pos)
    cno = ms2pos(i);

    for pno = cno   % (only exact scan here; neighbor logic kept commented)
        if pno<1 || pno>num_MS2, continue; end

        % If instrument varies across scans, adapt ActiveType/tol and K1 on the fly
        if length(unique(instruments))>1
            c_instrument = MS2_index(pno,6);
            if 3==c_instrument || 4==c_instrument
                ActiveType = 'ETD';
            else
                ActiveType = 'CID';
            end
            if mod(c_instrument,2)==1, tol = 0.4; else, tol = 0.02; end

            if 1==nhmass
                [K1,posn,posc] = get_key_ions1H(His,hno,Mods,ActiveType);
            else
                [K1,posn,posc] = get_key_ions1(His,hno,Mods,ActiveType);
            end
        end

        % Pull peaks for this MS2 scan
        IX    = index(pno):index(pno+1)-1;
        mz    = MS2_peaks(IX,1);
        inten = MS2_peaks(IX,2);

        % Basic noise filter: keep peaks >= 3× (scan’s noise estimate, col 8)
        I = find(inten >= 3*MS2_index(pno,8));
        mz = mz(I); inten = inten(I);

        % For each diagnostic ion, match nearest peak within tolerance
        for j=1:length(K1)
            ix1 = find(abs(mz-K1(j)) <= tol);
            [~,x1] = min(abs(mz(ix1)-K1(j))); %#ok<ASGLU>
            if ~isempty(ix1) && ms2intens(cno,j) < inten(ix1(x1))
                ms2intens(cno,j) = inten(ix1(x1));
            end
        end
    end
end
end
