function [cur_rts,cur_intens,cur_mono_isointens] = get_histone0(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,special)
%%
% GET_HISTONE0  Extract RT and integrated MS1 intensity for a histone peptide across charge states.
%
% Syntax:
%   [cur_rts, cur_intens, cur_mono_isointens] = get_histone0( ...
%       MS1_index, MS1_peaks, ptol, unitdiff, His, hno, special)
%
% Purpose (overview):
%   For a given histone peptide (row 'hno' in 'His'), this routine:
%     1) Builds MS1 extracted-ion profiles (XICs) around several isotopic m/z lines for the
%        FIRST charge state of the peptide (columns of His.pep_mz / His.pep_ch).
%     2) Determines a candidate apex scan/time near reference RTs with rules depending on
%        'special.ndebug' and on 'His.rt_ref' / (optionally) a file-derived RT for a known panel.
%     3) Integrates the monoisotopic chromatographic signal at that apex using GET_AREA
%        (which may deconvolve isotopic overlaps), yielding 'cur_intens(1)' and the apex 'cur_rts(1)'.
%     4) Repeats profile extraction and local integration for the REMAINING charge states, but
%        constrained to a narrow RT window centered at the first charge's apex (±delta).
%
% Inputs:
%   MS1_index : [nScans × ≥3] int/double; scan-wise index map into MS1_peaks.
%               Convention used here: MS1_index(k,2) = RT for scan k; MS1_index(k,3) = start row in MS1_peaks.
%   MS1_peaks : [nRows × 2] double; concatenated MS1 centroids as [m/z, intensity] across scans.
%   ptol      : scalar double; mass tolerance in ppm. (Special handling: ptol==100 → coerced to 10.)
%   unitdiff  : scalar double; isotopic mass spacing in Da (per ¹³C step), used as unitdiff/charge.
%   His       : struct with (at least) fields:
%                 • pep_mz  : [nPeptides × nCharge] target m/z per charge state
%                 • pep_ch  : [nPeptides × nCharge] charge states (z)
%                 • rt_ref  : [nPeptides × 1] reference RT per peptide
%                 • rt_unmod_orig : scalar or 1×1; reference for "unmodified" peptide RT
%                 • (optional) outfile : char; name of the current panel (used for special logic)
%                 • (optional) outpath : char; base path for auxiliary files (used in helper below)
%   hno       : scalar index selecting the peptide row in 'His'.
%   special   : struct with controls:
%                 • ndebug : 0 (normal), 1 (debug narrow window), -1 (normal but relaxed ref window).
%
% Outputs:
%   cur_rts            : [1 × nCharge] double; chosen apex RT per charge state (or fallback ref RT).
%   cur_intens         : [1 × nCharge] double; integrated area per charge (0 if not found).
%   cur_mono_isointens : [N × 1] double; the FULL monoisotopic XIC vector used for the FIRST charge state.
%
% Notes:
%   - This function does NOT alter 'His'; it only reads fields.
%   - Apex selection depends on 'special.ndebug' and on peptide-specific panels (e.g., 'H3_02_9_17').
%   - Area integration is delegated to GET_AREA, which performs overlap-aware integration if ptol<100.
%

[npep,ncharge] = size(His.pep_mz);%#ok
cur_rts = zeros([1,ncharge]);
cur_intens = zeros([1,ncharge]);
num_MS1 = size(MS1_index,1);
end_rt = MS1_index(num_MS1,2);
if end_rt>90
    gradient = 2;
else
    gradient = 1;
end;
% Heuristic gradient code used later by get_H33_27_40_unmod_rt to apply RT offsets.

% get MS1 profile
delta = 1;  % ±1 RT unit window for anchoring neighbors and for other charges.
c_mz = His.pep_mz(hno,1);
c_ch = His.pep_ch(hno,1);
c_ref_isomzs = [c_mz-unitdiff/c_ch c_mz c_mz+unitdiff/c_ch c_mz+2*unitdiff/c_ch];
if ptol==100
    ptol = 10;
end;
% Special convention: an input of 100 ppm is treated as 10 ppm (historical code behavior).

if ptol>100 && c_ch>=3
    nC13 = 1;
else
    nC13 = 0;
end;
% nC13 hints profile builder to consider 13C-related handling for high-charge states at coarse ppm.

if His.rt_unmod_orig==His.rt_ref(1) && 1~=special.ndebug
    rt_i1 = 1;
    rt_i2 = num_MS1;
else
    rt1 = His.rt_ref(hno)-5;
    rt2 = His.rt_ref(hno)+5;
    p = find( MS1_index(:,2)>=rt1 );
    rt_i1 = p(1);
    pp = find( MS1_index(:,2)<=rt2 );
    rt_i2 = pp(end);
end;
% Define scan-range [rt_i1:rt_i2] for initial search:
%   - default: whole run if 'unmod' reference equals His.rt_ref(1) and not in ndebug modes,
%   - otherwise: a ±5 min band around peptide-specific His.rt_ref(hno).

[c_isorts,c_ref_isointens] = GetProfiles(MS1_index,MS1_peaks,c_ref_isomzs,c_ch,ptol,nC13,rt_i1:rt_i2);
j = 2;
c_mono_isointens = c_ref_isointens(:,j);
cur_mono_isointens = c_mono_isointens;
% Keep the monoisotopic XIC vector for output (only for the first charge state).

% get rt and area
[nt,nb,top1_idx,inten_sum] = GetTopBottom(c_mono_isointens);
% GetTopBottom is assumed to:
%   - return candidate apex indices 'nt',
%   - segment boundaries 'nb',
%   - a primary index 'top1_idx',
%   - and a per-candidate intensity measure 'inten_sum'.

if 1==special.ndebug || -1==special.ndebug
    % Debug / relaxed-reference modes:
    if 1==special.ndebug
        llimit = -delta;
        rlimit = delta;
    else
        llimit = -2;
        if 1==isfield(His,'outfile') && 1==ismember(His.outfile,{'H3_02_9_17','H3_02a_9_17'})
            rlimit = 3;
        else
            rlimit = 2;
        end;
    end;
    ref_rt = His.rt_ref(hno);
    flag = zeros([1,length(nt)]);
    for i=1:length(nt)
        if llimit<=c_isorts(nt(i))-ref_rt && c_isorts(nt(i))-ref_rt<=rlimit
            flag(i) = 1;
        end;
    end;
    x = find(flag==1);
    if 0==isempty(x)
        %[tmp,id] = min(abs(c_isorts(nt(x))-ref_rt));%#ok
        if -1==special.ndebug && 1==isfield(His,'outfile') && 1==ismember(His.outfile,{'H3_02_9_17','H3_02a_9_17'}) && length(x)>=2
            [tmp,ix] = sort(inten_sum(x),'descend');%#ok
            ix = ix(1:2);
            [tmp,xx] = min( c_isorts(nt(x(ix))) );%#ok
            id = ix(xx);
        else
            [tmp,id] = max(inten_sum(x));%#ok
        end;
        top1_idx = x(id);
        cur_pos = nt(top1_idx);
        cur_rts(1) = c_isorts(cur_pos);
        cur_intens(1) = get_area(c_isorts,c_ref_isointens,nb,cur_pos,c_mz,c_ch,MS1_index,MS1_peaks,unitdiff,ptol);
    else
        cur_rts(1) = ref_rt;
        cur_intens(1) = 0;
    end;
else
    % Normal mode with special-case handling for specific panels:
    if His.rt_unmod_orig~=His.rt_ref(1) || ( 1==isfield(His,'outfile') && 1==strcmp(His.outfile,'H3_04v3_27_40') )
        if His.rt_unmod_orig~=His.rt_ref(1)
            ref_rt = His.rt_ref(1);
        elseif 1==isfield(His,'outfile') && 1==strcmp(His.outfile,'H3_04v3_27_40')
            ref_rt = get_H33_27_40_unmod_rt(His,gradient);
        end;
        flag = zeros([1,length(nt)]);
        for i=1:length(nt)
            if abs( c_isorts(nt(i))-ref_rt )<=delta
                flag(i) = 1;
            end;
        end;
        x = find(flag==1);
        if 0==isempty(x)
            %[tmp,id] = min(abs(c_isorts(nt(x))-ref_rt));%#ok
            [tmp,id] = max(inten_sum(x));%#ok
            top1_idx = x(id);
            cur_pos = nt(top1_idx);
            cur_rts(1) = c_isorts(cur_pos);
            cur_intens(1) = get_area(c_isorts,c_ref_isointens,nb,cur_pos,c_mz,c_ch,MS1_index,MS1_peaks,unitdiff,ptol);
        else
            cur_rts(1) = ref_rt;
            cur_intens(1) = 0;
        end;
    else
        if 0==isempty(nt) && c_isorts(nt(top1_idx))>4
            cur_pos = nt(top1_idx);
            cur_rts(1) = c_isorts(cur_pos);
            cur_intens(1) = get_area(c_isorts,c_ref_isointens,nb,cur_pos,c_mz,c_ch,MS1_index,MS1_peaks,unitdiff,ptol);
        else
            cur_rts(1) = His.rt_ref(hno);
            cur_intens(1) = 0;
        end;
    end;
end;
%{
if 0==cur_intens(1)
    return;
end;
%}
% Optional early return (disabled by design): if no intensity for the first charge, abort.

% Define a tight window ±delta around the chosen apex to process remaining charge states:
rt1 = cur_rts(1)-delta;
rt2 = cur_rts(1)+delta;
p = find( MS1_index(:,2)>=rt1 );
rt_i1 = p(1);
pp = find( MS1_index(:,2)<=rt2 );
rt_i2 = pp(end);

for jno=2:ncharge
    % get MS1 profile
    c_mz = His.pep_mz(hno,jno);
    c_ch = His.pep_ch(hno,jno);
    c_ref_isomzs = [c_mz-unitdiff/c_ch c_mz c_mz+unitdiff/c_ch c_mz+2*unitdiff/c_ch];
    if ptol>100 && c_ch>=3
        nC13 = 1;
    else
        nC13 = 0;
    end;
    [c_isorts,c_ref_isointens] = GetProfiles(MS1_index,MS1_peaks,c_ref_isomzs,c_ch,ptol,nC13,rt_i1:rt_i2);
    j = 2;
    c_mono_isointens = c_ref_isointens(:,j);

    % get rt and area
    [nt,nb,top1_idx,inten_sum] = GetTopBottom(c_mono_isointens);%#ok
    flag = zeros([1,length(nt)]);
    for i=1:length(nt)
        if abs( c_isorts(nt(i))-cur_rts(1) )<=delta
            flag(i) = 1;
        end;
    end;
    x = find(flag==1);
    if 0==isempty(x)
        %[tmp,id] = min(abs(c_isorts(nt(x))-cur_rts(1)));%#ok
        [tmp,id] = max(inten_sum(x));%#ok
        top1_idx = x(id);
        cur_pos = nt(top1_idx);
        cur_rts(jno) = c_isorts(cur_pos);
        cur_intens(jno) = get_area(c_isorts,c_ref_isointens,nb,cur_pos,c_mz,c_ch,MS1_index,MS1_peaks,unitdiff,ptol);
    else
        cur_rts(jno) = cur_rts(1);
        cur_intens(jno) = 0;
    end;
end;

function ref_rt = get_H33_27_40_unmod_rt(His,gradient)
%%
% GET_H33_27_40_UNMOD_RT  Panel-specific fallback RT for H3(27–40) unmodified, with gradient-based offset.
%
% Behavior:
%   - Reads '<His.outpath>/H3_04_27_40.xls', scans for a line '[rt]', then reads the
%     following lines to extract a numeric RT (tab-delimited).
%   - Applies an empirical offset based on 'gradient':
%       gradient==1 → ref_rt = rt0 - 0.3
%       gradient==2 → ref_rt = rt0 - 1.6
%   - If the file doesn't exist → fallback to His.rt_ref(1).

out_file0 = fullfile(His.outpath,'H3_04_27_40.xls');
if 0~=exist(out_file0,'file')
    fp = fopen(out_file0,'r');
    str = fgetl(fp);
    while 0==feof(fp) && 0==strcmp(str,'[rt]')
        str = fgetl(fp);
    end;
    str = fgetl(fp);%#ok
    str = fgetl(fp);
    fclose(fp);
    p = strfind(str,'	');
    rt0 = str2num( str(p(1)+1:p(2)-1) );%#ok
    if 1==gradient
        ref_rt = rt0-0.3;
    else
        ref_rt = rt0-1.6;
    end;
else
    ref_rt = His.rt_ref(1);
end;
