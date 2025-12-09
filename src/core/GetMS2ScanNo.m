function success = GetMS2ScanNo(ms2_fullfile)
%% GetMS2ScanNo — Parse a text MS2 export into compact MAT files (index + peaks)
%
% PLAIN-ENGLISH PURPOSE (for non-experts)
%   This function reads a vendor-neutral **text export** of MS2 spectra and creates:
%
%     <name>_MS2scans.mat  -> MS2_index (scan-wise metadata + pointers)
%     <name>_MS2peaks.mat  -> MS2_peaks (all m/z,intensity pairs concatenated)
%
%   It also relies on an **MS1 index** (created by GetMS1ScanNo) to backfill the
%   nearest MS1 scan/RT for each MS2 spectrum. If the MS2 MAT files already exist,
%   it returns immediately with success=1.
%
% EXPECTED INPUT FILE STRUCTURE (plain text, tab-separated headers)
%   Headers the parser looks for (literal TABs '\t' are part of the strings):
%     - 'H\tDataType'         ... then a char indicating data mode ('C' centroid / 'P' profile)
%     - 'S'                   ... new MS2 spectrum line; contains numeric fields (see below)
%     - 'I\tNumberOfPeaks'    ... (optional/original count; not used downstream)
%     - 'I\tRetTime'          ... MS2 retention time (usually seconds)
%     - 'I\tIonInjectionTime' ... (stored; not propagated in final table)
%     - 'I\tActivationType'   ... 'CID' | 'ETD' | 'HCD'
%     - 'I\tInstrumentType'   ... 'ITMS' | 'FTMS'
%     - 'I\tPrecursorScan'    ... MS1 scan number (link to MS1 table if present)
%     - 'I\tActivationCenter' ... isolation/activation center m/z (used to refine precursor m/z)
%     - 'I\tMonoiosotopicMz'  ... (misspelled in source; not used here)
%     - 'Z'                   ... charge/mass line; first occurrence parsed as: "<z>\t<MH+>"
%     - data lines: "<mz>\t<intensity>" repeated until next 'S' or EOF
%
% WHAT THE OUTPUT ARRAYS MEAN
%   MS2_index: [nScans × 8] double, columns are:
%       (1) MS1_scan     ... MS1 scan number associated to this MS2 (nearest/backfilled)
%       (2) ms2_rt       ... MS2 retention time (converted to minutes if values > 1000)
%       (3) MS2_scan     ... MS2 scan number as read from the 'S' line
%       (4) prec_mz      ... refined precursor m/z (aligned to ActivationCenter via isotopic offset)
%       (5) charge       ... precursor charge (from first 'Z' line)
%       (6) Fragtype     ... 1..6 encodes {CID-IT, CID-FT, ETD-IT, ETD-FT, HCD-IT, HCD-FT}
%       (7) ptr          ... *start pointer* into MS2_peaks for this MS2 scan
%       (8) baseline     ... per-scan baseline used to filter profile-like noise (0 if not applied)
%
%   MS2_peaks: [nPeaks × 2] double:
%       column 1 = m/z,  column 2 = intensity
%
%   To retrieve peaks for scan i:
%       start = MS2_index(i,7);
%       stop  = (i<nScans) ? MS2_index(i+1,7)-1 : size(MS2_peaks,1);
%       scan_peaks = MS2_peaks(start:stop, :);
%
% IMPORTANT BEHAVIORS
%   • The input must be **centroided**; if 'H\tDataType' says 'P' (profile), we abort.
%   • We correct the reported precursor m/z (from the 'S' line) by snapping it to the
%     closest **isotopic offset** that best matches 'ActivationCenter' within 20 ppm.
%     This helps align to the monoisotopic (or proper isotopologue) precursor.
%   • If duplicate MS2 scan numbers are detected back-to-back, the duplicate is skipped.
%
% RETURN
%   success = 1 on success (files created or found), 0 otherwise.

success = 0;

% ----- 0) If cached MAT files exist, return -----
[datapath,dataname] = fileparts(ms2_fullfile);
MS2_scanfile = fullfile(datapath,[dataname,'_MS2scans.mat']);
MS2_peakfile = fullfile(datapath,[dataname,'_MS2peaks.mat']);
if 0~=exist(MS2_scanfile,'file') && 0~=exist(MS2_peakfile,'file')
    success = 1;
    return;
end

% ----- 1) Sanity check input text file -----
if 0==exist(ms2_fullfile,'file')
    disp([ms2_fullfile,': does not exist!']);
    return;
end

fid=fopen(ms2_fullfile,'r');
if -1==fid
    disp(['can not open: ',ms2_fullfile]);
    return;
end

%% ----- 2) Preallocation & constants -----
pmass = 1.007276;          % proton mass used to compute MH+ from m/z and charge
maxpeaknum   = 1e4;        % max peaks per MS2 spectrum
maxMS2num    = 1.5e5;      % initial number of MS2 scans
totalpeaknum = 4e7;        % initial total peaks across all MS2 scans

% Final column layout (see header): [MS1_scan, ms2_rt, MS2_scan, prec_mz, z, Fragtype, ptr, baseline]
MS2_index = zeros([maxMS2num,8]);
MS2_peaks = zeros([totalpeaknum,2]);
fno = 0;                   % actual MS2 scan count
pkno = 0;                  % total stored peaks
oldms2scan = 0;            % to skip immediate duplicate MS2 scan numbers

% Load MS1 index (required to backfill MS1 scan/RT for each MS2)
% NOTE: expects sibling folder structure ".../<parent>/MS1/<dataname>_MS1scans.mat"
MS1_scanfile = fullfile(fileparts(datapath),'MS1',[dataname,'_MS1scans.mat']);
load(MS1_scanfile);        % provides MS1_index [scan_no, rt, ptr, baseline]
num_MS1 = size(MS1_index,1);

% Header keywords (literal TABs)
keyword0 = 'H	DataType';
keyword1 = 'S';
keyword2 = 'I	NumberOfPeaks';
keyword3 = 'I	RetTime';
keyword4 = 'I	IonInjectionTime';
keyword5 = 'I	ActivationType';
keyword6 = 'I	InstrumentType';
keyword7 = 'I	PrecursorScan';
keyword8 = 'I	ActivationCenter';
keyword9 = 'I	MonoiosotopicMz';  %#ok<NASGU> (kept for completeness; not used)
keyword10 = 'Z';
len0 = length(keyword0);
len1 = length(keyword1);
len2 = length(keyword2);   %#ok<NASGU>
len3 = length(keyword3);
len4 = length(keyword4);
len5 = length(keyword5);
len6 = length(keyword6);
len7 = length(keyword7);
len8 = length(keyword8);
len9 = length(keyword9);   %#ok<NASGU>
len10 = length(keyword10);

%% ----- 3) Detect data mode (must be centroid) -----
str=fgets(fid);
while feof(fid)==0 && 0==strcmp( str(1:len0),keyword0 )
    str=fgets(fid);
end
MS2_datamode = str(len0+2);
if 1==strcmp('P',MS2_datamode)
    fprintf(1,'MS2 is profile mode, convert to centroid mode first!\n');
    fclose(fid);
    return;
end

% Progress counter
ct_prt = 0;
fprintf(1,'MS2 scans: ');

%% ----- 4) Main parse loop -----
str = fgets(fid);
while 0==feof(fid)
    if 1==strcmp( str(1:len1),keyword1 )
        % --- 4.1 Parse the 'S' line numeric fields ---
        % Internal holding buffer (10 slots; see mapping below)
        tmp_datam = zeros(1,10);

        % The 'S' line carries: [MS2_scan, <possibly repeated MS2_scan>, cur_precursor_mz]
        % We keep: tmp_datam(2) = MS2_scan, tmp_datam(3) = current precursor m/z
        tmp_datam(1:3) = str2num(str(len1+2:end)); %#ok<ST2NM>

        % Skip immediate duplicate MS2 scan numbers (robustness against malformed exports)
        if 0==oldms2scan - tmp_datam(2)
            str=fgets(fid);
            while feof(fid)==0 && 0==strcmp( str(1:len1),keyword1 )
                str=fgets(fid);
            end
            continue;
        else
            oldms2scan = tmp_datam(2);
        end

        % Progress printing
        fno = fno + 1;
        fprintf(repmat('\b',[1,ct_prt]));
        ct_prt = fprintf('%i',fno);

        % --- 4.2 Optional/original #peaks (ignored) ---
        str=fgets(fid); %#ok<NASGU>

        % --- 4.3 Retention time (seconds) ---
        str=fgets(fid);
        if 1==strcmp( str(1:len3),keyword3 )
            tmp_datam(6) = str2num(str(len3+2:end)); %#ok<ST2NM>  % ms2_rt
        end

        % --- 4.4 IonInjectionTime (kept for completeness) ---
        str=fgets(fid);
        if 1==strcmp( str(1:len4),keyword4 )
            tmp_datam(8) = str2num(str(len4+2:end)); %#ok<ST2NM>
        end

        % --- 4.5 ActivationType -> numeric code 1..3 ---
        str=fgets(fid);
        if 1==strcmp( str(1:len5),keyword5 )
            tmp_str = str(len5+2:len5+4);
            if     1==strcmp( tmp_str,'CID' ), tmp_datam(9)  = 1;
            elseif 1==strcmp( tmp_str,'ETD' ), tmp_datam(9)  = 2;
            elseif 1==strcmp( tmp_str,'HCD' ), tmp_datam(9)  = 3;
            end
        end

        % --- 4.6 InstrumentType -> numeric code 1..2 ---
        str=fgets(fid);
        if 1==strcmp( str(1:len6),keyword6 )
            tmp_str = str(len6+2:len6+5);
            if     1==strcmp( tmp_str,'ITMS' ), tmp_datam(10) = 1;
            elseif 1==strcmp( tmp_str,'FTMS' ), tmp_datam(10) = 2;
            end
        end

        % Fragtype 1..6 (pairs activation × analyzer):
        %   1=CID-IT, 2=CID-FT, 3=ETD-IT, 4=ETD-FT, 5=HCD-IT, 6=HCD-FT
        tmp_datam(7) = (tmp_datam(9)-1)*2 + tmp_datam(10);

        % --- 4.7 PrecursorScan: link to nearest MS1 scan/RT ---
        str=fgets(fid);
        if 1==strcmp( str(1:len7),keyword7 )
            tmp_datam(1) = str2num(str(len7+2:end)); %#ok<ST2NM>  % MS1_scan
            % If missing, fallback to nearest MS1 scan prior to this MS2 scan number
            p = find(MS1_index(1:num_MS1,1) <= tmp_datam(2));
            if 1==isempty(p)
                fno = fno - 1;
                continue;
            end
            if 1==isempty(tmp_datam(1)) || 0==tmp_datam(1)
                tmp_datam(1) = MS1_index(p(end),1);
            end
            % Backfill MS2 RT from MS1 RT if needed
            tmp_datam(6) = MS1_index(p(end),2);
        end

        % --- 4.8 Activation center m/z (for refining precursor m/z) ---
        str=fgets(fid);
        if 1==strcmp( str(1:len8),keyword8 )
            cen_mz = str2num(str(len8+2:end)); %#ok<ST2NM>
        end

        % --- 4.9 MonoiosotopicMz (kept for compatibility; not used) ---
        str=fgets(fid); %#ok<NASGU>

        % --- 4.10 'Z' line(s): charge and MH+; compute refined precursor m/z ---
        str=fgets(fid);
        if 1==strcmp( str(1:len10),keyword10 )
            tmp_datam(4:5) = str2num(str(len10+2:end)); %#ok<ST2NM>   % [charge, MH+]
            % Snap precursor m/z to the isotopic offset that best matches ActivationCenter
            acc_mz = get_acc_mz(cen_mz, tmp_datam(3), tmp_datam(4));
            acc_mh = acc_mz*tmp_datam(4) - (tmp_datam(4)-1)*pmass;     % MH+ = m/z*z − (z−1)*H+
            tmp_datam(3) = acc_mz;   % refined precursor m/z
            tmp_datam(5) = acc_mh;   % MH+ for reference
        end
        % Skip any additional 'Z' lines (only first charge state is used)
        str=fgets(fid);
        while 1==strcmp( str(1:len10),keyword10 )
            tmp_datam(4:5) = [0 0];  % ignore subsequent charges
            str=fgets(fid);
        end

        % --- 4.11 Read centroided peaks until next 'S' or EOF ---
        mz    = zeros([1,maxpeaknum]);
        inten = zeros([1,maxpeaknum]);
        pnum  = 0;
        while feof(fid)==0 && 0==strcmp( str(1:len1),keyword1 )
            if ~('0'<=str(1) && str(1)<='9')
                str = fgets(fid);
                continue;
            end
            pnum = pnum + 1;
            tmp = sscanf(str,'%f',2);
            mz(pnum)    = tmp(1);
            inten(pnum) = tmp(2);
            str = fgets(fid);
        end
        if 1==feof(fid)
            % include last line at EOF if it contains a pair
            pnum = pnum + 1;
            tmp = sscanf(str,'%f',2);
            mz(pnum)    = tmp(1);
            inten(pnum) = tmp(2);
        end

        % Drop non-positive intensities
        IX = find(inten>0);
        mz    = mz(IX);
        inten = inten(IX);

        % --- 4.12 (Optional) centroiding — not implemented (expects centroid input) ---

        % --- 4.13 Baseline filter (if enough points); discard degenerate spectra ---
        if 1==isempty(mz) || mz(end)-mz(1)<10
            fno = fno - 1;
            continue;
        end
        if length(inten)>100
            baseline = GetBaseline(inten);
            IX = find(inten>=baseline);
            mz    = mz(IX);
            inten = inten(IX);
        else
            baseline = 0;
        end
        npk = length(mz);

        % --- 4.14 Store metadata and append peaks ---
        % Reorder fields into final layout: [MS1_scan, ms2_rt, MS2_scan, prec_mz, z, Fragtype]
        x = [1 6 2 3 4 7];
        new_datam = tmp_datam(x);
        MS2_index(fno,1:8) = [new_datam npk baseline];
        MS2_peaks(pkno+1:pkno+npk,1) = mz;
        MS2_peaks(pkno+1:pkno+npk,2) = inten;
        pkno = pkno + npk;

    else
        % Not an 'S' line -> continue reading
        str = fgets(fid);
    end
end
fclose(fid);
fprintf(repmat('\b',[1,ct_prt]));
fprintf('%i',fno);
fprintf(1,'\n');

%% ----- 5) Finalize & save -----
% Trim preallocations
if fno<maxMS2num
    IX = 1:fno;
    MS2_index = MS2_index(IX,:);
end
if pkno<totalpeaknum
    IX = 1:pkno;
    MS2_peaks = MS2_peaks(IX,:); %#ok<NASGU>
end

% Convert column 7 from "npk per scan" to *cumulative start pointers* (+1)
tmp = MS2_index(1:fno,7);
MS2_index(1:fno,7) = cumsum(tmp) + 1;

% Convert MS2 RT to minutes if values look like seconds (heuristic threshold)
if MS2_index(fno,2)>1000
    MS2_index(1:fno,2) = MS2_index(1:fno,2)/60; %#ok<NASGU>
end

% Save to MAT files
save(MS2_scanfile,'MS2_index');
save(MS2_peakfile,'MS2_peaks');

success = 1;


% ====================== nested helpers ======================
function baseline = GetBaseline(inten)
%% GetBaseline — Mode-of-log-intensity baseline estimator (per MS2 spectrum)
%   Convert intensities to log10, histogram them with coarse bins, take the
%   mode bin center as baseline, and return 10^(mode center).
loginten = log10(inten);
t = min(loginten):0.08:max(loginten);
[n,xout] = hist(loginten,t);
[~,idx] = max(n);
baseline = 10^xout(idx);


function acc_mz = get_acc_mz(cen_mz,cur_mz,cur_chg)
%% get_acc_mz — Snap precursor m/z to an isotopic offset near ActivationCenter
%   Given:
%     cen_mz  ... 'ActivationCenter' m/z (from header)
%     cur_mz  ... precursor m/z reported on the 'S' line
%     cur_chg ... precursor charge z (from 'Z' line)
%
%   We generate candidate precursor m/z values by shifting cur_mz by ±k * (1.0032 / z)
%   (13C spacing scaled by charge) for k in [-1..8], keep those within ±20 ppm of
%   cen_mz, and select the nearest. If none within tolerance, return cen_mz.
ptol = 20;                 % ppm tolerance
unitdiff = 1.0032;         % ~13C − 12C mass difference
sets = [-1 0 1 2 3 4 5 6 7 8];
mzs = cur_mz + sets*unitdiff/cur_chg;
ix = find(abs(mzs-cen_mz)<ptol*cen_mz*1e-6);
if ~isempty(ix)
    [~,i] = min(abs( mzs(ix)-cen_mz ));
    acc_mz = mzs(ix(i));
else
    acc_mz = cen_mz;
end
