function GetPSM2(cur_outpath,out_filename,His,pep_rts,pep_intens,isorts,mono_isointens,MS1_index,MS1_peaks,MS2_index,ptol,unitdiff)
%% GetPSM2 — Alternate PSM export (labeling focused on mod_short)
%
% PLAIN-ENGLISH PURPOSE (for non-experts)
%   Same high-level goal as GetPSM: gather peptide-spectrum matches (PSMs)
%   by intersecting MS1-derived retention-time windows with MS2 precursor
%   evidence (m/z & charge). The **difference** is how labels/fields are
%   composed:
%     • psm.pep     now uses:  "<mod_short_stripped>_<pep_seq>"
%     • psm.seq     now stores: "<mod_short_stripped>"   (not the bare pep_seq)
%     • psm.mod0    now stores: "<pep_seq>"             (swapped role vs GetPSM)
%   This is useful when you want the visual ID to prioritize the annotated
%   modified sequence (mod_short) rather than the raw peptide sequence.
%
% INPUTS (same shapes as in GetPSM; see that header for details)
%   cur_outpath, out_filename
%   His:   .pep_seq (char), .pep_mz (npep×nchg), .pep_ch (npep×nchg),
%          .mod_short{npep×1}, .mod_type{npep×1}
%   pep_rts (npep×1), pep_intens (npep×1)
%   isorts (nMS1×1), mono_isointens (nMS1×npep)
%   MS1_index, MS1_peaks (from GetMS1ScanNo)
%   MS2_index (from GetMS2ScanNo)
%   ptol (ppm; if 100 -> uses 10), unitdiff (~1.0032 Da)
%
% OUTPUT FILES (written into <cur_outpath>/psm/)
%   <out_filename>.mat     — MATLAB struct 'psm' (alternate labeling)
%   <out_filename>.plabel  — pLabel-compatible text file
%
% KEY BEHAVIOR DIFFERENCES vs GetPSM
%   • Labeling:
%       - We derive c_seq := His.mod_short{*} stripped after last '.'
%         (e.g., "H3K9ac.ASAPSTK" -> "ASAPSTK"). If no '.' exists, take it as-is.
%       - psm.pep{ino} = "<c_seq>_<His.pep_seq>"
%       - PSM rows: psm.seq = c_seq; psm.mod0 = His.pep_seq
%     Rationale: make the "modified peptide" string primary in displays.
%
% MATCHING LOGIC (unchanged)
%   • Define RT window per peptide around MS1 apex (GetTopBottom + GetLocal).
%   • Inside that RT range, keep MS2 whose (charge matches) and
%     (precursor m/z within ±ptol ppm of theoretical), allowing M, M+1, M+2
%     by adding ±k*unitdiff/z to the theoretical m/z.
%   • If measured precursor m/z is out-of-ppm, snap to the strongest MS1
%     peak within the ppm window on the linked MS1 scan (cosmetic correction).
%
% SAFETY NOTES
%   • This function expects consistent His.mod_short entries. If they are
%     not in the "prefix.string" format, stripping by last '.' will simply
%     keep the whole string (no error, but labeling may differ from intent).
%   • Tolerance quirk: ptol==100 is treated as 10 ppm (historical behavior).

    % ----- 0) Tolerance normalization & output folder -----
    if ptol==100
        ptol = 10;
    end
    psm_outpath = fullfile(cur_outpath,'psm');
    if 0==exist(psm_outpath,'dir') && 0==mkdir(psm_outpath)
        fprintf(1,'can not create: %s\n',psm_outpath);
        return;
    end
    out_file1 = fullfile(psm_outpath,[out_filename,'.mat']);
    out_file2 = fullfile(psm_outpath,[out_filename,'.plabel']);

    % ----- 1) Peptide subset (RT sanity check) -----
    [npep,nchg] = size(His.pep_mz);
    ix = find(pep_rts(1:npep,1)>4);   % ignore very-early eluters (<4 min)
    if isempty(ix)
        return;
    end

    % ----- 2) RT windows ("terminus") from monoisotopic XIC -----
    % Same approach as GetPSM: locate MS1 scan near apex, then refine a
    % left/right scan range around local maximum if peptide has signal.
    nplot    = length(ix);
    terminus = repmat(0,[nplot,2]);
    for ino=1:nplot
        cno = ix(ino);
        p = find(isorts<=pep_rts(cno,1));
        c_ms1pos = p(end);                        % MS1 scan index at apex (from the left)
        c_mono_isointens = mono_isointens(:,cno); % MS1 monoisotopic trace

        if pep_intens(cno,1)>0
            [nt,nb] = GetTopBottom(c_mono_isointens); %#ok (segment)
            [localmax_rt,localmax_inten,IX] = GetLocal(c_ms1pos,isorts,c_mono_isointens,nb); %#ok
            if isempty(IX)
                terminus(ino,1:2) = [c_ms1pos c_ms1pos];
            else
                terminus(ino,1:2) = [IX(1) IX(end)];
            end
        else
            terminus(ino,1:2) = [c_ms1pos c_ms1pos];
        end
    end

    % ----- 3) Paths and labels -----
    [layoutpath,cur_raw] = fileparts(fileparts(cur_outpath));
    datapath  = fileparts(layoutpath);
    ms2_path  = fullfile(datapath,'MS2');
    p0 = strfind(cur_raw,'_');
    cur_rawname = cur_raw(p0(1)+1:end);  % remove prefix up to first underscore
    p = strfind(out_filename,'_');
    c_prot = out_filename(1:p(1)-1);     % protein/histone label prefix
    sets = [0 1 2];                      % allow M, M+1, M+2
    num_MS2 = size(MS2_index,1);
    num_MS1 = size(MS1_index,1);
    index = [1;MS1_index(1:num_MS1,3)];  % MS1 pointers

    % ----- 4) Initialize PSM container (alternate labeling) -----
    % Build psm.pep = "<mod_short_stripped>_<pep_seq>"
    psm.pep = repmat({''},[npep,1]);
    for ino=1:npep
        c_seq = His.mod_short{ino};
        p = strfind(c_seq,'.');             % we keep piece after the last '.'
        if ~isempty(p)
            c_seq = c_seq(p(end)+1:end);
        end
        psm.pep{ino,1} = [c_seq,'_',His.pep_seq];
    end
    psm.nms2 = repmat(0,[npep,1]);

    % flat PSM lists (rows per accepted MS2)
    psm.fname = {};
    psm.prot  = {};
    psm.seq   = {};   % <-- stores c_seq (mod_short stripped)
    psm.mod0  = {};   % <-- stores His.pep_seq (raw sequence)
    psm.mod1  = {};
    psm.mod2  = {};
    psm.emz0  = [];
    psm.emz   = [];
    psm.tmz   = [];
    psm.chg   = [];
    psm.rt    = [];
    fno = 0;

    % ----- 5) Collect MS2 matches per RT window -----
    for ino=1:nplot
        if terminus(ino,1)==terminus(ino,2)
            continue; % degenerate window -> skip
        end
        rt1 = isorts(terminus(ino,1));
        rt2 = isorts(terminus(ino,2));

        % MS2 within RT window
        p  = find( MS2_index(:,2)>=rt1 );
        if isempty(p), continue; end
        i1 = p(1);
        pp = find( MS2_index(:,2)<=rt2 );
        if isempty(pp), continue; end
        i2 = pp(end);

        % Check precursor match: same charge AND within ppm (allowing ±k*unitdiff/z)
        cno  = ix(ino);
        flag = repmat(0,[num_MS2,1]);   % flag(i)=column jno in pep_mz/pep_ch that matched
        for i=i1:i2
            cen_mz = MS2_index(i,4);
            cen_ch = MS2_index(i,5);
            for jno=1:nchg
                c_mz = His.pep_mz(cno,jno);
                c_ch = His.pep_ch(cno,jno);
                if c_ch~=cen_ch, continue; end
                mzs = c_mz + sets*unitdiff/c_ch;
                IX = find(abs(mzs-cen_mz)<ptol*cen_mz*1e-6); %#ok
                if ~isempty(IX)
                    flag(i) = jno;  % remember which charge column matched
                    break;
                end
            end
        end

        % Precompute the c_seq (mod_short stripped) for this peptide index
        c_seq = His.mod_short{cno};
        p = strfind(c_seq,'.');
        if ~isempty(p)
            c_seq = c_seq(p(end)+1:end);
        end

        % ----- 6) Emit PSM rows if there are matches -----
        ms2pos = find(flag>0);
        if ~isempty(ms2pos)              % note: code uses "if 0==isempty(ms2pos)" (same effect)
            psm.nms2(cno,1) = length(ms2pos);

            for no=1:length(ms2pos)
                i = ms2pos(no);
                cur_scan = MS2_index(i,3);
                cur_chg  = MS2_index(i,5);

                fno = fno + 1;
                psm.fname{fno,1} = [cur_rawname,'.',num2str(cur_scan),'.',num2str(cur_scan),'.',num2str(cur_chg)];
                psm.prot{fno,1}  = c_prot;

                % *** Main differences vs GetPSM ***
                psm.seq{fno,1}  = c_seq;              % primary display = mod_short stripped
                psm.mod0{fno,1} = His.pep_seq;        % store raw sequence here
                psm.mod1{fno,1} = His.mod_type{cno};  % long mod string (text)
                psm.mod2{fno,1} = change_modtype(c_seq,His.mod_type{cno}); % numeric mapping for pLabel

                % m/z fields (with possible MS1-based correction if out-of-ppm)
                psm.emz0(fno,1) = MS2_index(i,4);
                psm.emz(fno,1)  = MS2_index(i,4);
                c_tmz = His.pep_mz(cno,flag(i));
                psm.tmz(fno,1)  = c_tmz;

                if abs(psm.emz(fno,1)-c_tmz)>ptol*c_tmz*1e-6
                    % Find MS1 scan at or before this MS2 scan number
                    p = find( MS1_index(1:num_MS1,1)<=cur_scan );
                    cur_ms1pos = p(end);

                    % Read that MS1 spectrum
                    IX    = index(cur_ms1pos):index(cur_ms1pos+1)-1;
                    mz    = MS1_peaks(IX,1);
                    inten = MS1_peaks(IX,2);

                    % Snap to strongest MS1 peak in ppm window around c_tmz
                    cur_ptol = ptol*c_tmz*1e-6;
                    left  = c_tmz - cur_ptol;
                    right = c_tmz + cur_ptol;
                    x = find( mz>=left & mz<=right ); %#ok
                    if ~isempty(x)
                        [~,xx] = max(inten(x));
                        psm.emz(fno,1) = mz(x(xx));
                    end
                end

                psm.chg(fno,1) = cur_chg;
                psm.rt(fno,1)  = pep_rts(cno,1);
            end
        end
    end

    % ----- 7) Save outputs -----
    save(out_file1,'psm');
    output_plabel(out_file2,ms2_path,cur_rawname,psm);
end


% ====================== helpers ======================

function mod2 = change_modtype(pep,mod1)
%% change_modtype — Map "pos,mod;pos,mod;..." into pLabel numeric IDs
% Same logic as in GetPSM. Note: here 'pep' receives c_seq (mod_short stripped).
%
% INPUT
%   pep  : display peptide (here, the mod_short-stripped sequence)
%   mod1 : "pos,code;pos,code;..." (pos=0 => peptide N-terminus)
% OUTPUT
%   mod2 : "pos,ID;pos,ID;..." where ID enumerates:
%           1 = pr[PEP_N] (reserved), then GetMods() over names × residue sets.

    mod2 = '';

    pos1 = strfind(mod1,',');
    if isempty(pos1)
        return;  % no pairs to translate
    end

    Mods = GetMods();
    pos2 = [0 strfind(mod1,';')];
    for jno = 1 : length(pos1)
        cpos = str2num( mod1(pos2(jno)+1:pos1(jno)-1) ); %#ok<ST2NM>
        cmod = mod1(pos1(jno)+1:pos2(jno+1)-1);

        if 0==cpos
            % Special-case: position 0 + 'pr' => ID 1 ('pr[PEP_N]')
            if strcmp(cmod,'pr')  % ++++
              mod2 = [mod2,' 0,1']; %#ok<AGROW>
            end
        else
            % Walk the catalog to find numeric ID
            bflag = 0;
            no = 1;
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
%% output_plabel — Write a minimal .plabel for pLabel (same structure as GetPSM)
% The only difference is that spectrum blocks will reflect the alternate
% labeling (psm.seq is 'c_seq' here). Mod catalog remains identical.

    fp2 = fopen(out_file2,'w');
    if -1==fp2
        disp(['can not open the file: ',out_file2]);
        return;
    end

    % File header
    fprintf(fp2,'[FilePath]\r\n');
    fprintf(fp2,'File_Path=%s\r\n',fullfile(ms2_path,[cur_rawname,'.ms2']));

    % Modification catalog
    Mods = GetMods();
    no = 1;
    fprintf(fp2,'[Modification]\r\n');
    fprintf(fp2,'1=pr[PEP_N]\r\n'); % ++++ reserved entry
    for i=1:length(Mods.name)
        for j=1:length(Mods.set{i})
            no = no + 1;
            fprintf(fp2,'%d=%s[%s]\r\n',no,lower(Mods.name{i}),Mods.set{i}(j));
        end
    end

    % No cross-links
    fprintf(fp2,'[xlink]\r\n');
    fprintf(fp2,'xlink=NULL\r\n');

    % Totals
    fprintf(fp2,'[Total]\r\n');
    fprintf(fp2,'total=%d\r\n',sum(psm.nms2));

    % Spectrum blocks
    for ino = 1 : sum(psm.nms2)
        fprintf(fp2,'[Spectrum%d]\r\n',ino);
        fname = upper(psm.fname{ino,1});
        fprintf(fp2,'name=%s\r\n',fname);
        fprintf(fp2,'pep1=0 %s 1 %s\r\n',psm.seq{ino,1},psm.mod2{ino,1});
    end

    fclose(fp2);
end
