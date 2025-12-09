function output_histone(cur_outpath,out_filename,His,pep_intens,pep_rts)
%%
% OUTPUT_HISTONE
% Purpose:
%   Export per-peptide, per-charge intensities and retention times to a
%   tab-separated .xls (human/Excel-friendly), and also persist a .mat
%   containing a compact AUC-like table (RT, Area, GlobalFraction) for
%   the subset of peptides marked as display==1 in His.display.
%
% Inputs:
%   cur_outpath : output directory (string)
%   out_filename: base filename (string, no extension)
%   His         : structure describing the peptide set and annotation with:
%                 - pep_seq   : char, peptide sequence (header)
%                 - mod_short : cellstr (npep×1), short labels per peptide
%                 - mod_type  : cellstr (npep×1), detailed mod spec
%                 - pep_mz    : double (npep×ncharge), m/z per charge
%                 - pep_ch    : double/int (npep×ncharge), charge values
%                 - rt_ref    : double (npep×1), reference RTs
%                 - display   : logical/double (npep×1), 1=visible/export to .mat
%   pep_intens  : double (npep×ncharge), per-peptide per-charge intensities
%   pep_rts     : double (npep×ncharge), per-peptide per-charge RTs
%
% Outputs:
%   Creates two files under cur_outpath:
%     1) <out_filename>.xls  (tab-separated text with [area] and [rt] blocks)
%     2) <out_filename>.mat  (variables: auc [npep2×3], His [possibly reduced])
%
% Notes:
%   - The .xls is tab-separated text, not a binary Excel XLS. Excel/LibreOffice
%     open it fine; TSV is used for portability.
%   - Fractions are safeguarded with eps to avoid division-by-zero warnings.

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

% Open file handle for writing (text mode)
fp = fopen(out_file1,'w');
if -1==fp
    % If file cannot be opened, report and abort (no .mat will be written)
    fprintf(1,'can not open: %s\n',out_file1);
    return;
end

% Header: peptide sequence (single line)
fprintf(fp,'%s\r\n',His.pep_seq);

% ---- [area] block ----
fprintf(fp,'[area]\r\n');
fprintf(fp,'peptide\t');
for jno=1:ncharge
    % Two columns per charge: absolute area and within-row fraction
    fprintf(fp,'(+%d)\t\t',His.pep_ch(1,jno));
end
fprintf(fp,'total\t\r\n');

% Body: per peptide row, per charge area and fractions, then row total and global fraction
for ino=1:npep
    fprintf(fp,'%s\t',His.mod_short{ino});
    for jno=1:ncharge
        % Absolute area in scientific notation; within-row fraction as float
        fprintf(fp,'%e\t%f\t',pep_intens(ino,jno),pep_intens(ino,jno)/(eps+row_intens(ino)));
    end
    % Row total and global fraction (relative to sum of all row totals)
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

% RT per peptide and charge (two decimals)
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
% Write the .mat (compact AUC)
% ---------------------------

% Select only peptides flagged for display==1
ix = find(His.display==1);
if 1==isempty(ix)
    % Nothing to export to .mat (but the .xls was written); exit early.
    return;
end
npep2 = length(ix);

% auc: [RT(first charge), SumArea(all charges), GlobalFraction]
auc = zeros([npep2,3]);
for ino=1:npep2
    % RT: take the first charge's RT (column 1)
    auc(ino,1) = pep_rts(ix(ino),1);
    % Area: sum across all available charges
    auc(ino,2) = sum(pep_intens(ix(ino),1:ncharge));
end
total_intens = sum(auc(1:npep2,2));
for ino=1:npep2
    % Global fraction normalized over the display==1 subset
    auc(ino,3) = auc(ino,2)/(eps+total_intens);
end

% If not all peptides are displayed, reduce His to the visible subset for consistency
if npep2<npep
    His = reduce_His(His,ix); %#ok
end

% Persist variables
save(out_file2,'auc','His');

end % function output_histone


function His0 = reduce_His(His,ix)
%%
% REDUCE_HIS
% Purpose:
%   Create a consistent subset of the His structure restricted to rows in ix.
%   All dependent fields are subset accordingly, preserving shapes.
%
% Inputs:
%   His : full His struct
%   ix  : indices to keep (column vector or row vector)
%
% Output:
%   His0: reduced His struct with the same field semantics

npep2 = length(ix);

% Keep the peptide sequence as-is (single sequence header)
His0.pep_seq = His.pep_seq;

% mod_short (cellstr)
for ino=1:npep2
    His0.mod_short{ino,1} = His.mod_short{ix(ino),1};
end

% mod_type (cellstr)
for ino=1:npep2
    His0.mod_type{ino,1} = His.mod_type{ix(ino),1};
end

% pep_mz and pep_ch (double matrices); respect the number of charge columns
for ino=1:npep2
    His0.pep_mz(ino,1:size(His.pep_ch,2)) = His.pep_mz(ix(ino),1:size(His.pep_ch,2));
    His0.pep_ch(ino,1:size(His.pep_ch,2)) = His.pep_ch(ix(ino),1:size(His.pep_ch,2));
end

% rt_ref (column vector)
for ino=1:npep2
    His0.rt_ref(ino,1) = His.rt_ref(ix(ino),1);
end

end % function reduce_His
