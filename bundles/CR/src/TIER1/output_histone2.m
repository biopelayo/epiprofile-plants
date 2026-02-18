function output_histone2(cur_outpath,out_filename,His,pep_intens,pep_rts,onlyme)
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
auc = zeros([npep,3]);
for ino=1:npep
    auc(ino,1) = pep_rts(ino,1);
    auc(ino,2) = sum(pep_intens(ino,1:ncharge));
end;
total_intens = sum(auc(1:npep,2));
for ino=1:npep
    auc(ino,3) = auc(ino,2)/(eps+total_intens);
end;
[His,auc] = combineLH(cur_outpath,out_filename,His,auc,onlyme);%#ok
save(out_file2,'auc','His');

function [His1,auc1] = combineLH(cur_outpath,out_filename,His2,auc2,onlyme)
%%

[path1,name1] = fileparts(cur_outpath);
[path2,name2] = fileparts(path1);
out_file1 = fullfile(fullfile(fileparts(path2),name2,name1),[out_filename,'.mat']);
if 0==exist(out_file1,'file')
    fprintf(1,'%s: not exist.\n',out_file1);
    His1 = His2;
    auc1 = auc2;
    return;
end;
load(out_file1);
npep1 = length(His.mod_type);
if 1==onlyme;
    flag = zeros(npep1,1);
    for ino=1:npep1
        if 1==strcmp(out_filename,'H3_08_117_128') || 1==strcmp(out_filename,'H4_06_79_92')
            pos = strfind(His.mod_short{ino},'unmod');
        else
            pos = strfind(His.mod_type{ino},'me');
        end;
        if 0==isempty(pos)
            flag(ino) = 1;
        end;
    end;
    ix = find(flag==1);
    npep1 = length(ix);
    His1.pep_seq = His.pep_seq;
    His1.mod_short = His.mod_short(ix);
    His1.mod_type = His.mod_type(ix);
    His1.pep_mz(1:npep1,1:size(His.pep_ch,2)) = His.pep_mz(ix,1:size(His.pep_ch,2));
    His1.pep_ch(1:npep1,1:size(His.pep_ch,2)) = His.pep_ch(ix,1:size(His.pep_ch,2));
    His1.rt_ref = His.rt_ref(ix);
    auc1 = auc(ix,1:3);
else
    His1 = His;
    auc1 = auc;
end;

npep2 = length(His2.mod_type);
for ino=1:npep2
    His1.mod_short{ino+npep1,1} = His2.mod_short{ino,1};
    His1.mod_type{ino+npep1,1} = His2.mod_type{ino,1};
    His1.pep_mz(ino+npep1,1:size(His2.pep_ch,2)) = His2.pep_mz(ino,1:size(His2.pep_ch,2));
    His1.pep_ch(ino+npep1,1:size(His2.pep_ch,2)) = His2.pep_ch(ino,1:size(His2.pep_ch,2));
    His1.rt_ref(ino+npep1,1) = His2.rt_ref(ino,1);
    auc1(ino+npep1,1:3) = auc2(ino,1:3);
end;
npep = npep1+npep2;
total_intens = sum(auc1(1:npep,2));
for ino=1:npep
    auc1(ino,3) = auc1(ino,2)/(eps+total_intens);
end;