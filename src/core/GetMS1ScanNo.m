function success = GetMS1ScanNo(ms1_fullfile)
%% GetMS1ScanNo — Parse a text MS1 export into compact MAT files (index + peaks)
%
% PLAIN-ENGLISH PURPOSE (non-experts)
%   This function reads a vendor-independent **text export** of MS1 scans and
%   produces two compact MATLAB files used throughout the pipeline:
%
%     <name>_MS1scans.mat  -> MS1_index (scan-wise metadata + pointers), MS1Type
%     <name>_MS1peaks.mat  -> MS1_peaks (all m/z,intensity pairs concatenated)
%
%   It also handles caching: if those MAT files already exist (and contain
%   MS1Type), it will return immediately with success=1.
%
% EXPECTED INPUT FILE FORMAT (plain text)
%   The parser looks for simple *line headers* and data lines:
%     - 'H\tDataType'        ... then a single char indicating data mode ('C'/'P')
%     - 'S'                  ... start of a new MS1 scan (followed by scan number)
%     - 'I\tRetTime'         ... retention time (RT), usually in seconds
%     - 'I\tInstrumentType'  ... instrument code (e.g., 'FTMS', 'ITMS')
%     - data lines: "<mz><TAB><intensity>" repeated until the next 'S' or EOF
%
% WHAT THE OUTPUT ARRAYS MEAN
%   MS1_index: [nScans × 4] double, columns are:
%       (1) scan_no      ... scan identifier as read from the file
%       (2) rt           ... retention time (converted to minutes if >1000)
%       (3) ptr          ... *start pointer* into MS1_peaks for this scan
%       (4) baseline     ... per-scan baseline (0 if not computed)
%   MS1_peaks: [nPeaks × 2] double:
%       column 1 = m/z,  column 2 = intensity
%
%   For scan i, its peaks live in:
%       start = MS1_index(i,3);
%       stop  = (i<nScans) ? MS1_index(i+1,3)-1 : size(MS1_peaks,1);
%       scan_peaks = MS1_peaks(start:stop, :);
%
% SAFETY / ROBUSTNESS
%   - If the input is *profile mode* (DataType='P'), we abort and ask to centroid
%     first (this code expects centroided MS1).
%   - A per-scan baseline is estimated (histogram of log10 intensities) when the
%     scan has >100 points; peaks below baseline are discarded.
%
% RETURN VALUE
%   success = 1 on success (MAT files created or already present), 0 otherwise.

    success = 0;

    % ----- 0) Check for cached MAT files -----
    [datapath,dataname] = fileparts(ms1_fullfile);
    MS1_scanfile = fullfile(datapath,[dataname,'_MS1scans.mat']);
    MS1_peakfile = fullfile(datapath,[dataname,'_MS1peaks.mat']);
    if 0~=exist(MS1_scanfile,'file') && 0~=exist(MS1_peakfile,'file')
        % For legacy files that may miss MS1Type, add a default 'FTMS' once
        load(MS1_scanfile);
        if 0==exist('MS1Type') %#ok
            MS1Type = 'FTMS'; %#ok<NASGU>
            save(MS1_scanfile,'MS1_index','MS1Type');
        end
        success = 1;
        return;
    end

    % ----- 1) Check input text file -----
    if 0==exist(ms1_fullfile,'file')
        disp([ms1_fullfile,': does not exist!']);
        return;
    end
    fid = fopen(ms1_fullfile,'r');
    if -1==fid
        disp([ms1_fullfile,': can not open!']);
        return;
    end

    % ----- 2) Preallocation & constants -----
    % Heuristics for big files; adjust if needed for very large runs
    maxpeaknum    = 1e4;   % max peaks per *single* MS1 scan (m/z,int pairs)
    maxMS1num     = 1.5e5; % initial number of MS1 scans
    totalpeaknum  = 4e7;   % initial total peaks across all scans

    % MS1_index columns: [scan_no, rt, ptr/npk (temp), baseline]
    MS1_index = zeros([maxMS1num,4]);
    MS1_peaks = zeros([totalpeaknum,2]);  % [m/z, intensity]
    fno = 0;   % actual number of scans encountered
    pkno = 0;  % running count of stored peaks

    % Keywords (note the literal TAB '\t' in strings)
    keyword0 = 'H	DataType';
    keyword1 = 'S';
    keyword2 = 'I	RetTime';
    keyword3 = 'I	InstrumentType';
    len0 = length(keyword0);
    len1 = length(keyword1);
    len2 = length(keyword2);
    len3 = length(keyword3);

    % ----- 3) Read header to get data mode ('C' centroid / 'P' profile) -----
    str=fgets(fid);
    while feof(fid)==0 && 0==strcmp( str(1:len0),keyword0 )
        str=fgets(fid);
    end
    MS1_datamode = str(len0+2);  % read the char after "H<TAB>DataType<TAB>"
    if 1==strcmp('P',MS1_datamode)
        fprintf(1,'MS1 is profile mode, convert to centroid mode first!\n');
        fclose(fid);
        return;
    end

    % Progress indicator (prints a running scan count)
    ct_prt = 0;
    fprintf(1,'MS1 scans: ');

    % ----- 4) Main parse loop: scan-by-scan -----
    str = fgets(fid);
    while 0==feof(fid)
        if 1==strcmp( str(1:len1),keyword1 )
            % New scan begins
            fno = fno + 1;
            fprintf(repmat('\b',[1,ct_prt]));   % backspace previous counter
            ct_prt = fprintf('%i',fno);

            % 4.1 Scan number (integer after 'S')
            scan_no = str2num(str(len1+2:end)); %#ok<ST2NM>
            scan_no = scan_no(1);

            % 4.2 Retention time (line starting with 'I<tab>RetTime')
            str = fgets(fid);
            while feof(fid)==0 && 0==strcmp( str(1:len2),keyword2 )
                str=fgets(fid);
            end
            rt_no = str2double(str(len2+2:end));  % RT as provided (usually seconds)

            % 4.3 IonInjectionTime (skipped)
            str = fgets(fid); %#ok<NASGU>

            % 4.4 InstrumentType (read on first scan and saved as 4-char code)
            str = fgets(fid);
            if 1==fno
                MS1Type = str(len3+2:len3+5); %#ok<NASGU>
                % If you wanted to abort on low-res ITMS, you'd check and return here.
            end

            % 4.5 Read peaks until next 'S' or EOF
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

            % First, drop non-positive intensities
            IX = find(inten>0);
            mz    = mz(IX);
            inten = inten(IX);

            % 4.6 (optional) Centroiding — not implemented here (expects centroid input)

            % 4.7 Baseline & filtering (skip if scan is too small)
            if 1==isempty(mz) || mz(end)-mz(1)<10
                % Likely malformed scan; discard and roll back scan counter
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

            % 4.8 Store scan metadata & append peaks
            % Temporarily put npk in column 3; we will convert to pointers later
            MS1_index(fno,1:4) = [scan_no rt_no npk baseline];
            MS1_peaks(pkno+1:pkno+npk,1) = mz;
            MS1_peaks(pkno+1:pkno+npk,2) = inten;
            pkno = pkno + npk;
        else
            % Not a scan header -> read next line
            str = fgets(fid);
        end
    end
    fclose(fid);
    fprintf(repmat('\b',[1,ct_prt]));
    fprintf('%i',fno);
    fprintf(1,'\n');

    % ----- 5) Finalize arrays -----
    % Trim preallocations
    if fno<maxMS1num
        IX = 1:fno;
        MS1_index = MS1_index(IX,:);
    end
    if pkno<totalpeaknum
        IX = 1:pkno;
        MS1_peaks = MS1_peaks(IX,:); %#ok<NASGU>
    end

    % Convert column 3 from "npk per scan" to *cumulative start pointers* (+1)
    tmp = MS1_index(1:fno,3);
    MS1_index(1:fno,3) = cumsum(tmp) + 1;

    % Convert RT to minutes if values look like seconds (heuristic threshold)
    if MS1_index(fno,2)>1000 % e.g., > ~16.7 min in seconds
        MS1_index(1:fno,2) = MS1_index(1:fno,2)/60; %#ok<NASGU>
    end

    % ----- 6) Save to MAT files (cache) -----
    save(MS1_scanfile,'MS1_index','MS1Type');
    save(MS1_peakfile,'MS1_peaks');

    success = 1;
end


% ====================== nested helpers ======================
function baseline = GetBaseline(inten)
%% GetBaseline — Mode-of-log-intensity baseline estimator (per scan)
%   Rationale:
%     - Convert intensities to log10 and histogram them.
%     - The *mode* bin (highest count) corresponds to the noise floor / baseline.
%     - Return 10^(mode bin center) as the baseline threshold.
    loginten = log10(inten);
    t = min(loginten):0.08:max(loginten);   % coarse bins; 0.08 ~ 10^(0.08)=1.2 fold
    [n,xout] = hist(loginten,t);
    [~,idx] = max(n);
    baseline = 10^xout(idx);
end
