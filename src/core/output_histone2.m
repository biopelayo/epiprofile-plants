function output_histone2(cur_outpath,out_filename,His,pep_intens,pep_rts,onlyme)
%%
% OUTPUT_HISTONE2
% Purpose:
%   1) Export a tab-separated ".xls" (readable by Excel/LibreOffice) with
%      per-charge areas and intra-row/global fractions, and per-charge RTs.
%   2) Build an AUC-like table [RT(first charge), SumArea(all charges), GlobalFraction].
%   3) Optionally COMBINE the current His/auc with a pre-existing .mat that
%      lives under a "sibling" directory inferred from cur_outpath (see combineLH).
%   4) Re-normalize global fractions AFTER combination so the sum of col-2 equals 1
%      (modulo eps) across the combined set.
%
% Inputs:
%   cur_outpath : output directory for the current run
%   out_filename: base file name (no extension)
%   His         : struct describing peptides (fields: pep_seq, mod_short, mod_type,
%                 pep_mz, pep_ch, rt_ref; no display filter applied here)
%   pep_intens  : double [npep×ncharge], areas/intensities per peptide and charge
%   pep_rts     : double [npep×ncharge], RT per peptide and charge
%   onlyme      : scalar flag:
%                 - if 1, combineLH will subset the LOADED reference set to methylations
%                   (substring 'me' in mod_type), except for 2 special filenames where
%                   it will subset 'unmod' (substring in mod_short).
%
% Outputs (files written):
%   <cur_outpath>/<out_filename>.xls  (TSV-like text with [area] and [rt])
%   <cur_outpath>/<out_filename>.mat  (variables: auc [n×3], His [combined])
%
% Notes:
%   - The .xls is plain text with tabs; it's not a binary XLS.
%   - RT in auc(:,1) uses the FIRST charge column as-is.

out_file1 = fullfile(cur_outpath,[out_filename,'.xls']);
out_file2 = fullfile(cur_outpath,[out_filename,'.mat']);
[npep,ncharge] = size(His.pep_mz);

% ---------------------------
% Write the .xls (tabulated)
% ---------------------------

% Row-wise totals across all charge states
row_intens = zeros([npep,1]);
for ino=1:npep
    row_intens(ino) = sum(pep_intens(ino,:));
end
total_intens = sum(row_intens);

% Open handle
fp = fopen(out_file1,'w');
if -1==fp
    % If cannot open, report and abort (no .mat will be written)
    fprintf(1,'can not open: %s\n',out_file1);
    return;
end

% Header: peptide sequence (single line)
fprintf(fp,'%s\r\n',His.pep_seq);

% ---- [area] block ----
fprintf(fp,'[area]\r\n');
fprintf(fp,'peptide\t');
for jno=1:ncharge
    % Two columns per charge: absolute area, within-row fraction
    fprintf(fp,'(+%d)\t\t',His.pep_ch(1,jno));
end
fprintf(fp,'total\t\r\n');

% Body: per peptide row
for ino=1:npep
    fprintf(fp,'%s\t',His.mod_short{ino});
    for jno=1:ncharge
        % Absolute area in scientific notation; fraction within the peptide's row
        fprintf(fp,'%e\t%f\t',pep_intens(ino,jno),pep_intens(ino,jno)/(eps+row_intens(ino)));
    end
    % Row total and global fraction (relative to total_intens)
    fprintf(fp,'%e\t%f\r\n',row_intens(ino),row_intens(ino)/(eps+total_intens));
end

% ---- [rt] block ----
fprintf(fp,'\r\n[rt]\r\n');
fprintf(fp,'peptide\t');
for jno=1:ncharge
    if jno==ncharge
        fprintf(fp,'(+%d)',His.pep_ch(1,jno));
    else
        fprintf(fp,'(+%d)\t',His.pep_ch(1,jno));
    end
end
fprintf(fp,'\r\n');

% RT table (two decimals)
for ino=1:npep
    fprintf(fp,'%s\t',His.mod_short{ino});
    for jno=1:ncharge
        if jno==ncharge
            fprintf(fp,'%.2f',pep_rts(ino,jno));
        else
            fprintf(fp,'%.2f\t',pep_rts(ino,jno));
        end
    end
    fprintf(fp,'\r\n');
end

% Close the tabulated file
fclose(fp);

% ---------------------------
% Build auc and COMBINE
% ---------------------------

% AUC-like table: [ RT(first charge), SumArea(all charges), GlobalFraction ]
auc = zeros([npep,3]);
for ino=1:npep
    % RT: first charge column (as provided upstream)
    auc(ino,1) = pep_rts(ino,1);
    % Area: sum across all charges
    auc(ino,2) = sum(pep_intens(ino,1:ncharge));
end
% Global fractions over the CURRENT set
total_intens = sum(auc(1:npep,2));
for ino=1:npep
    auc(ino,3) = auc(ino,2)/(eps+total_intens);
end

% Combine with a sibling reference .mat if present (and possibly subset it)
[His,auc] = combineLH(cur_outpath,out_filename,His,auc,onlyme); %#ok

% Persist the combined result
save(out_file2,'auc','His');

end % function output_histone2



function [His1,auc1] = combineLH(cur_outpath,out_filename,His2,auc2,onlyme)
%%
% COMBINELH
% Purpose:
%   Locate a "reference" .mat under a sibling directory (derived by going up
%   two levels from cur_outpath and reconstructing a parallel path), load it,
%   optionally SUBSET it (onlyme), then APPEND the current set (His2/auc2)
%   and re-normalize the global fractions.
%
% Inputs:
%   cur_outpath : current output directory
%   out_filename: base name to look for in the reference directory
%   His2, auc2  : current set to be appended after the (optional) reference
%   onlyme      : if 1, subset the LOADED set based on:
%                 - if out_filename is 'H3_08_117_128' or 'H4_06_79_92':
%                      keep rows with 'unmod' in mod_short
%                 - else:
%                      keep rows with 'me' in mod_type (methylations)
%
% Outputs:
%   His1, auc1  : merged (reference[maybe subset] + current), with auc1(:,3)
%                 re-normalized to the total sum of column 2.

[path1,name1] = fileparts(cur_outpath);           % name1 ~ current histone folder
[path2,name2] = fileparts(path1);                 % name2 ~ parent folder name
out_file1 = fullfile(fullfile(fileparts(path2),name2,name1),[out_filename,'.mat']);
% out_file1 points to: parent_of_parent(path2)/name2/name1/out_filename.mat
% i.e., a "sibling" tree mirroring <name1>/<out_filename>.mat at a parallel root.

if 0==exist(out_file1,'file')
    % If no reference file, report and fall back to current set
    fprintf(1,'%s: not exist.\n',out_file1);
    His1 = His2;
    auc1 = auc2;
    return;
end

% Load reference 'His' and 'auc' (these variable names are expected inside the .mat)
load(out_file1);
npep1 = length(His.mod_type);

% Optionally subset the LOADED set (reference) before merge
if 1==onlyme;
    flag = zeros(npep1,1);
    for ino=1:npep1
        % Two special filenames: keep 'unmod' by checking mod_short
        if 1==strcmp(out_filename,'H3_08_117_128') || 1==strcmp(out_filename,'H4_06_79_92')
            pos = strfind(His.mod_short{ino},'unmod');
        else
            % General case: keep methylated forms (substring 'me' in mod_type)
            pos = strfind(His.mod_type{ino},'me');
        end
        if 0==isempty(pos)
            flag(ino) = 1;
        end
    end
    ix = find(flag==1);
    npep1 = length(ix);

    % Rebuild His1 and auc1 from the filtered reference set
    His1.pep_seq = His.pep_seq;
    His1.mod_short = His.mod_short(ix);
    His1.mod_type  = His.mod_type(ix);
    His1.pep_mz(1:npep1,1:size(His.pep_ch,2)) = His.pep_mz(ix,1:size(His.pep_ch,2));
    His1.pep_ch(1:npep1,1:size(His.pep_ch,2)) = His.pep_ch(ix,1:size(His.pep_ch,2));
    His1.rt_ref = His.rt_ref(ix);
    auc1 = auc(ix,1:3);
else
    % No subsetting: take the loaded reference as-is
    His1 = His;
    auc1 = auc;
end

% Append the CURRENT set (His2/auc2) after the (possibly filtered) reference set
npep2 = length(His2.mod_type);
for ino=1:npep2
    His1.mod_short{ino+npep1,1} = His2.mod_short{ino,1};
    His1.mod_type{ino+npep1,1}  = His2.mod_type{ino,1};
    His1.pep_mz(ino+npep1,1:size(His2.pep_ch,2)) = His2.pep_mz(ino,1:size(His2.pep_ch,2));
    His1.pep_ch(ino+npep1,1:size(His2.pep_ch,2)) = His2.pep_ch(ino,1:size(His2.pep_ch,2));
    His1.rt_ref(ino+npep1,1)    = His2.rt_ref(ino,1);
    auc1(ino+npep1,1:3)         = auc2(ino,1:3);
end

% Recompute global fractions over the MERGED set
npep = npep1+npep2;
total_intens = sum(auc1(1:npep,2));
for ino=1:npep
    auc1(ino,3) = auc1(ino,2)/(eps+total_intens);
end

end % function combineLH
