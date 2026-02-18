function output_histone(cur_outpath,out_filename,His,pep_intens,pep_rts)
%%

out_file1 = fullfile(cur_outpath,[out_filename,'.xls']);
out_file2 = fullfile(cur_outpath,[out_filename,'.mat']);
[npep,ncharge] = size(His.pep_mz);

% output out_file1
% write intens
row_intens = zeros([npep,1]);
for ino=1:npep
    row_intens(ino) = sum(pep_intens(ino,:));
end;
total_intens = sum(row_intens);

fp = fopen(out_file1,'w');
if -1==fp
    fprintf(1,'can not open: %s\n',out_file1);
    return;
end;
fprintf(fp,'%s\r\n',His.pep_seq);
fprintf(fp,'[area]\r\n');
fprintf(fp,'peptide\t');
for jno=1:ncharge
    fprintf(fp,'(+%d)\t\t',His.pep_ch(1,jno));
end;
fprintf(fp,'total\t\r\n');

for ino=1:npep
    fprintf(fp,'%s\t',His.mod_short{ino});
    for jno=1:ncharge
        fprintf(fp,'%e\t%f\t',pep_intens(ino,jno),pep_intens(ino,jno)/(eps+row_intens(ino)));
    end;
    fprintf(fp,'%e\t%f\r\n',row_intens(ino),row_intens(ino)/(eps+total_intens));
end;

% write rt
fprintf(fp,'\r\n[rt]\r\n');
fprintf(fp,'peptide\t');
for jno=1:ncharge
    if jno==ncharge
        fprintf(fp,'(+%d)',His.pep_ch(1,jno));
    else
        fprintf(fp,'(+%d)\t',His.pep_ch(1,jno));
    end;
end;
fprintf(fp,'\r\n');

for ino=1:npep
    fprintf(fp,'%s\t',His.mod_short{ino});
    for jno=1:ncharge
        if jno==ncharge
            fprintf(fp,'%.2f',pep_rts(ino,jno));
        else
            fprintf(fp,'%.2f\t',pep_rts(ino,jno));
        end;
    end;
    fprintf(fp,'\r\n');
end;

fclose(fp);

% output out_file2
% RT Area Ratio
ix = find(His.display==1);
if 1==isempty(ix)
    return;
end;
npep2 = length(ix);
auc = zeros([npep2,3]);
for ino=1:npep2
    auc(ino,1) = pep_rts(ix(ino),1);
    auc(ino,2) = sum(pep_intens(ix(ino),1:ncharge));
end;
total_intens = sum(auc(1:npep2,2));
for ino=1:npep2
    auc(ino,3) = auc(ino,2)/(eps+total_intens);
end;
if npep2<npep
    His = reduce_His(His,ix);%#ok
end;
save(out_file2,'auc','His');

function His0 = reduce_His(His,ix)
%%

npep2 = length(ix);
His0.pep_seq = His.pep_seq;
for ino=1:npep2
    His0.mod_short{ino,1} = His.mod_short{ix(ino),1};
end;
for ino=1:npep2
    His0.mod_type{ino,1} = His.mod_type{ix(ino),1};
end;
for ino=1:npep2
    His0.pep_mz(ino,1:size(His.pep_ch,2)) = His.pep_mz(ix(ino),1:size(His.pep_ch,2));
    His0.pep_ch(ino,1:size(His.pep_ch,2)) = His.pep_ch(ix(ino),1:size(His.pep_ch,2));
end;
for ino=1:npep2
    His0.rt_ref(ino,1) = His.rt_ref(ix(ino),1);
end;