function GetBenchmark(cur_outpath)
%%

psm_outpath = fullfile(cur_outpath,'psm');

benchmarkfile = fullfile(psm_outpath,'identification_list.xls');
fp = fopen(benchmarkfile,'w');
if -1==fp
    fprintf('can not open:%s\n',benchmarkfile);
    return;
end;

fprintf(fp,'filename\tmeasured m/z\tcalculated m/z\tcharge\tppm\tsequence\tmodification1\tmodification2\thistone type\tretention time(min)\n');
matfiles = dir(fullfile(psm_outpath,'*.mat'));
for ino=1:length(matfiles)
    cur_file = fullfile(psm_outpath,matfiles(ino).name);
    load(cur_file);
    nlen = length(psm.fname);
    for jno=1:nlen
        ppm = 1e6*(psm.emz(jno)-psm.tmz(jno))/psm.tmz(jno);
        fprintf(fp,'%s\t%.6f\t%.6f\t%d\t%.6f\t%s\t%s\t%s\t%s\t%.4f\n',psm.fname{jno},psm.emz(jno),psm.tmz(jno),psm.chg(jno),ppm,psm.seq{jno},psm.mod0{jno},psm.mod1{jno},psm.prot{jno},psm.rt(jno));
    end;
end;
fclose(fp);