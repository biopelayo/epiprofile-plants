function GetPSM(cur_outpath,out_filename,His,pep_rts,pep_intens,isorts,mono_isointens,MS1_index,MS1_peaks,MS2_index,ptol,unitdiff)
%% GetPSM — Build a minimal PSM table and a pLabel .plabel file from MS1+MS2 evidence
%
% PLAIN-ENGLISH PURPOSE (for non-experts)
%   This function collects peptide-spectrum matches (PSMs) for a single protein/
%   peptide sequence (contained in "His") by:
%     1) defining where each peptide elutes in RT (a start–end "terminus" window)
%        using the monoisotopic MS1 chromatogram,
%     2) scanning MS2 spectra that fall inside that RT window,
%     3) keeping only MS2 whose precursor m/z&charge match the peptide within ppm,
%     4) exporting a compact MATLAB struct 'psm' (*.mat) and a text *.plabel file
%        compatible with pLabel/pFind-style tools.
%
% HIGH-LEVEL INPUTS (typical shapes; adapt to your data model)
%   cur_outpath     : output folder (this function will create ./psm/ inside)
%   out_filename    : base name for outputs (without extension)
%   His             : peptide/protein descriptor struct with (per-peptide arrays):
%       .pep_seq        (char)   -> peptide sequence (same for all charge states)
%       .pep_mz (npep×nchg)      -> theoretical precursor m/z per charge
%       .pep_ch (npep×nchg)      -> integer charge states aligned with pep_mz
%       .mod_short{npep×1}       -> short mod string per peptide (e.g., "K9ac;...") for labels
%       .mod_type {npep×1}       -> long mod string "pos,mod;pos,mod;..." (used for pLabel)
%   pep_rts         : [npep×1] apex/representative RT (min) per peptide
%   pep_intens      : [npep×1] summary intensity per peptide (used to decide windows)
%   isorts          : [nMS1×1] RT (min) per MS1 scan
%   mono_isointens  : [nMS1×npep] monoisotopic XIC per peptide (MS1 traces)
%   MS1_index       : [nMS1×4]  columns: [scan_no, rt, ptr, baseline]  (from GetMS1ScanNo)
%   MS1_peaks       : [nPeaks×2] concatenated [m/z, intensity]          (from GetMS1ScanNo)
%   MS2_index       : [nMS2×8]  columns: [MS1_scan, rt, MS2_scan, prec_mz, z, Fragtype, ptr, baseline]
%   ptol            : ppm tolerance for precursor matching (special case: 100 -> use 10)
%   unitdiff        : isotopic spacing for 13C (~1.0032 Da) used to allow M, M+1, M+2 checks
%
% OUTPUTS (files written into <cur_outpath>/psm/)
%   <out_filename>.mat     -> MATLAB struct 'psm' with minimal PSM fields for downstream use
%   <out_filename>.plabel  -> text file listing spectra and peptide+mods for pLabel visualization
%
% KEY IDEAS
%   • RT window ("terminus") per peptide is derived from local MS1 behavior around its apex:
%     we call GetTopBottom (to segment) and GetLocal (to refine) on the monoisotopic XIC.
%   • MS2 hits are accepted only if (precursor_charge matches) AND (precursor m/z within ppm
%     of theoretical m/z, allowing ±k*(unitdiff/z) for isotopic picking).
%   • If the MS2 precursor m/z is slightly off the theoretical window, we "snap" it to the
%     closest MS1 peak within tolerance on the *linked* MS1 scan (cosmetic fix for reporting).
%   • Mod encoding: we convert EpiProfile-style "pos,mod;..." into pLabel numeric codes
%     using GetMods() to enumerate the mod catalog (see helper 'change_modtype' below).

    % ----- 0) Tolerance normalization & output folder -----
    if ptol==100
        ptol = 10;  % historical knob: 100 means "use 10 ppm"
    end
    psm_outpath = fullfile(cur_outpath,'psm');
    if 0==exist(psm_outpath,'dir') && 0==mkdir(psm_outpath)
        fprintf(1,'can not create: %s\n',psm_outpath);
        return;
    end
    out_file1 = fullfile(psm_outpath,[out_filename,'.mat']);
    out_file2 = fullfile(psm_outpath,[out_filename,'.plabel']);

    % ----- 1) Peptide subset to process (pep RT must be sensible > 4 min) -----
    [npep,nchg] = size(His.pep_mz);
    ix = find(pep_rts(1:npep,1)>4);
    if isempty(ix)
        % nothing to report
        return;
    end

    % ----- 2) Define RT windows ("terminus") per peptide via local MS1 profile -----
    % For each peptide with valid RT:
    %   • find MS1 scan index of its apex RT,
    %   • if peptide is "present" (pep_intens>0), refine a left/right boundary around
    %     that apex using GetTopBottom (coarse valleys) + GetLocal (trim zeros; ensure
    %     peak is included). Otherwise, set a degenerate window [apex, apex].
    nplot = length(ix);
    terminus = zeros([nplot,2]);   % columns are [left_scan_index, right_scan_index]
    for ino=1:nplot
        cno = ix(ino);
        p = find(isorts<=pep_rts(cno,1));
        c_ms1pos = p(end);                          % MS1 scan index nearest RT apex (from the left)
        c_mono_isointens = mono_isointens(:,cno);   % monoisotopic XIC across scans

        if pep_intens(cno,1)>0
            [nt,nb] = GetTopBottom(c_mono_isointens);            %#ok (segment trace into peaks)
            [localmax_rt,localmax_inten,IX] = GetLocal(c_ms1pos,isorts,c_mono_isointens,nb); %#ok
            if isempty(IX)
                terminus(ino,1:2) = [c_ms1pos c_ms1pos];
            else
                terminus(ino,1:2) = [IX(1) IX(end)];             % use refined scan-range around local max
            end
        else
            terminus(ino,1:2) = [c_ms1pos c_ms1pos];
        end
    end

    % ----- 3) Resolve paths to MS2 raw (for .plabel header) and protein label -----
    [layoutpath,cur_raw] = fileparts(fileparts(cur_outpath));  % assume: <...>/<RAW>/<psm>
    datapath = fileparts(layoutpath);
    ms2_path = fullfile(datapath,'MS2');
    p0 = strfind(cur_raw,'_');
    cur_rawname = cur_raw(p0(1)+1:end);           % strip prefix up to first underscore
    p = strfind(out_filename,'_');
    c_prot = out_filename(1:p(1)-1);              % protein/histone code to display
    sets = [0 1 2];                                % allow M, M+1, M+2 matching

    num_MS2 = size(MS2_index,1);
    num_MS1 = size(MS1_index,1);
    index = [1;MS1_index(1:num_MS1,3)];            % MS1 pointers to peaks

    % ----- 4) Initialize PSM container (MAT struct) -----
    % psm.pep: per-peptide label "<sequence>_<shortmod>"
    psm.pep = repmat({''},[npep,1]);
    for ino=1:npep
        psm.pep{ino,1} = [His.pep_seq,'_',His.mod_short{ino}];
    end
    psm.nms2 = zeros([npep,1]);   % how many MS2 match each peptide (count)

    % flat PSM lists (one row per accepted MS2)
    psm.fname = {};   % spectrum label "<raw>.<ms2scan>.<ms2scan>.<z>"
    psm.prot  = {};   % protein/histone name (e.g., "H3")
    psm.seq   = {};   % peptide sequence (no mods in-sequence)
    psm.mod0  = {};   % short mod string (as in His.mod_short{pep})
    psm.mod1  = {};   % long mod string (as in His.mod_type{pep})
    psm.mod2  = {};   % numeric pLabel mapping "pos,ID;pos,ID;..." (built below)
    psm.emz0  = [];   % measured precursor m/z (original from MS2_index)
    psm.emz   = [];   % (possibly corrected) precursor m/z (snapped to MS1 if out of ppm)
    psm.tmz   = [];   % theoretical precursor m/z used for matching
    psm.chg   = [];   % charge
    psm.rt    = [];   % peptide RT (min; from pep_rts)

    fno = 0;          % running PSM row index

    % ----- 5) For each peptide RT window, scan MS2 and keep matching precursors -----
    for ino=1:nplot
        % skip degenerate windows (no width)
        if terminus(ino,1)==terminus(ino,2)
            continue;
        end
        rt1 = isorts(terminus(ino,1));
        rt2 = isorts(terminus(ino,2));

        % MS2 within the peptide RT window
        p = find( MS2_index(:,2)>=rt1 );
        if isempty(p), continue; end
        i1 = p(1);
        pp = find( MS2_index(:,2)<=rt2 );
        if isempty(pp), continue; end
        i2 = pp(end);

        % Check precursor m/z±isotope and charge for each MS2 in window
        cno = ix(ino);                      % peptide index in His arrays
        flag = zeros([num_MS2,1]);          % 0=no match; j>0 -> matched charge-state column j in pep_mz/pep_ch
        for i=i1:i2
            cen_mz = MS2_index(i,4);        % measured precursor m/z
            cen_ch = MS2_index(i,5);        % measured charge

            for jno=1:nchg
                c_mz = His.pep_mz(cno,jno); % theoretical precursor m/z for this peptide & charge
                c_ch = His.pep_ch(cno,jno);
                if c_ch~=cen_ch
                    continue;               % charge mismatch -> skip
                end

                % Allow isotopic offsets: M (+0), M+1, M+2 (scaled by z)
                mzs = c_mz + sets*unitdiff/c_ch;
                IX = find(abs(mzs-cen_mz)<ptol*cen_mz*1e-6); %#ok
                if ~isempty(IX)
                    flag(i) = jno;          % remember which charge-state column matched
                    break;
                end
            end
        end

        % ----- 6) Materialize PSM rows for all matched MS2 in the RT window -----
        ms2pos = find(flag>0);
        if ~isempty(ms2pos)
            psm.nms2(cno,1) = length(ms2pos);

            for no=1:length(ms2pos)
                i = ms2pos(no);
                cur_scan = MS2_index(i,3);    % MS2 scan number
                cur_chg  = MS2_index(i,5);    % charge (z)

                fno = fno + 1;
                % psm.fname format: "<raw>.<ms2scan>.<ms2scan>.<z>"
                psm.fname{fno,1} = [cur_rawname,'.',num2str(cur_scan),'.',num2str(cur_scan),'.',num2str(cur_chg)];
                psm.prot{fno,1}  = c_prot;
                psm.seq{fno,1}   = His.pep_seq;
                psm.mod0{fno,1}  = His.mod_short{cno};
                psm.mod1{fno,1}  = His.mod_type{cno};
                psm.mod2{fno,1}  = change_modtype(His.pep_seq,His.mod_type{cno}); % numeric pLabel encoding

                psm.emz0(fno,1)  = MS2_index(i,4);                 % reported precursor m/z
                psm.emz(fno,1)   = MS2_index(i,4);                 % (may be corrected below)

                c_tmz = His.pep_mz(cno,flag(i));                   % theoretical m/z (the matched charge col)
                psm.tmz(fno,1)   = c_tmz;

                % If measured m/z is outside ppm tolerance, snap to nearest MS1 peak
                if abs(psm.emz(fno,1)-c_tmz) > ptol*c_tmz*1e-6
                    % find the MS1 scan *at or before* this MS2 scan number
                    p = find( MS1_index(1:num_MS1,1)<=cur_scan );
                    cur_ms1pos = p(end);

                    % Extract that MS1 spectrum
                    IX = index(cur_ms1pos):index(cur_ms1pos+1)-1;
                    mz    = MS1_peaks(IX,1);
                    inten = MS1_peaks(IX,2);

                    % Window around theoretical m/z
                    cur_ptol = ptol*c_tmz*1e-6;
                    left  = c_tmz - cur_ptol;
                    right = c_tmz + cur_ptol;
                    x = find( mz>=left & mz<=right );

                    % Replace emz with the m/z of the strongest MS1 peak in-window
                    if ~isempty(x)
                        [~,xx] = max(inten(x)); %#ok
                        psm.emz(fno,1) = mz(x(xx));
                    end
                end

                psm.chg(fno,1) = cur_chg;
                psm.rt(fno,1)  = pep_rts(cno,1);
            end
        end
    end

    % ----- 7) Save PSM struct and emit .plabel file -----
    save(out_file1,'psm');
    output_plabel(out_file2,ms2_path,cur_rawname,psm);
end


% ====================== helpers ======================

function mod2 = change_modtype(pep,mod1)
%% change_modtype — Map "pos,mod;pos,mod;..." into pLabel numeric IDs
% INPUT
%   pep  : peptide sequence (string)
%   mod1 : long mod string "pos,code;pos,code;..." where pos=0 means N-term mod
%
% OUTPUT
%   mod2 : pLabel-style "pos,ID;pos,ID;..." where ID enumerates the mod catalog:
%          ID=1 is reserved for 'pr[PEP_N]' (N-term propionyl, here extended as 0,1),
%          subsequent IDs follow GetMods() order across names×sets.
%
% LOGIC
%   • If pos=0 and mod='pr', emit "0,1" (N-term propionyl special-case).
%   • Otherwise, find (name,set) in the mod catalog (GetMods), where 'set' restricts
%     the allowed residue letters (e.g., 'K', 'STY', ...). The numeric ID is 1-based:
%       1 -> N-term pr[PEP_N], then iterate mods across all residues to assign IDs.

    mod2 = '';

    pos1 = strfind(mod1,',');
    if isempty(pos1)
        return;  % no "pos,mod" pairs
    end

    Mods = GetMods();
    pos2 = [0 strfind(mod1,';')];
    for jno = 1 : length(pos1)
        cpos = str2num( mod1(pos2(jno)+1:pos1(jno)-1) ); %#ok<ST2NM>
        cmod = mod1(pos1(jno)+1:pos2(jno+1)-1);

        if 0==cpos
            % Special-case: position 0 (peptide N-term in this code) and 'pr'
            if strcmp(cmod,'pr')     % ++++ (historical special handling)
                mod2 = [mod2,' 0,1']; %#ok<AGROW>
            end
        else
            % Find numeric ID 'no' by walking the mod catalog in pLabel order
            bflag = 0;
            no = 1;  % 1 is 'pr[PEP_N]'; then we enumerate all residue-scoped mods
            for i=1:length(Mods.name)
                for j=1:length(Mods.set{i})
                    no = no + 1;
                    if strcmp(cmod,Mods.name{i}) && pep(cpos)==Mods.set{i}(j)
                        bflag = 1; break;
                    end
                end
                if bflag, break; end
            end
            mod2 = [mod2,' ',num2str(cpos),',',num2str(no)]; %#ok<AGROW>
        end
    end

    mod2 = mod2(2:end);   % strip leading space
end


function output_plabel(out_file2,ms2_path,cur_rawname,psm)
%% output_plabel — Write a minimal .plabel file for pLabel visualization
% FORMAT OVERVIEW
%   [FilePath]
%   File_Path=<full path to .ms2>
%
%   [Modification]
%   1=pr[PEP_N]
%   2=ac[K]
%   3=ac[S]
%   ...
%
%   [xlink]
%   xlink=NULL
%
%   [Total]
%   total=<number of spectra>
%
%   [Spectrum1]
%   name=<RAW.MS2SCAN.MS2SCAN.Z>   (uppercased)
%   pep1=0 <PEPSEQ> 1 <mod2>
%
% NOTES
%   • 'mod2' is the numeric mapping created in 'change_modtype'.
%   • 'File_Path' points to "<ms2_path>/<cur_rawname>.ms2" (adjust folder layout if needed).

    fp2 = fopen(out_file2,'w');
    if -1==fp2
        disp(['can not open the file: ',out_file2]);
        return;
    end

    % File path to the .ms2 (raw text spectra)
    fprintf(fp2,'[FilePath]\r\n');
    fprintf(fp2,'File_Path=%s\r\n',fullfile(ms2_path,[cur_rawname,'.ms2']));

    % Modifications catalog: ID=1 reserved for 'pr[PEP_N]', then enumerate GetMods
    Mods = GetMods();
    no = 1;
    fprintf(fp2,'[Modification]\r\n');
    fprintf(fp2,'1=pr[PEP_N]\r\n'); % ++++ special entry for N-term propionyl
    for i=1:length(Mods.name)
        for j=1:length(Mods.set{i})
            no = no + 1;
            fprintf(fp2,'%d=%s[%s]\r\n',no,lower(Mods.name{i}),Mods.set{i}(j));
        end
    end

    % No cross-links in this workflow
    fprintf(fp2,'[xlink]\r\n');
    fprintf(fp2,'xlink=NULL\r\n');

    % Total spectra reported
    fprintf(fp2,'[Total]\r\n');
    fprintf(fp2,'total=%d\r\n',sum(psm.nms2));

    % One block per spectrum; 'name' uppercased, 'pep1' holds sequence+mods
    for ino = 1 : sum(psm.nms2)
        fprintf(fp2,'[Spectrum%d]\r\n',ino);
        fname = upper(psm.fname{ino,1});
        fprintf(fp2,'name=%s\r\n',fname);
        fprintf(fp2,'pep1=0 %s 1 %s\r\n',psm.seq{ino,1},psm.mod2{ino,1});
    end

    fclose(fp2);
end
