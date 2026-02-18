function GetPSM(cur_outpath,out_filename,His,pep_rts,pep_intens,isorts,mono_isointens,MS1_index,MS1_peaks,MS2_index,ptol,unitdiff)
%%

if ptol==100
    ptol = 10;
end;
psm_outpath = fullfile(cur_outpath,'psm');
if 0==exist(psm_outpath,'dir') && 0==mkdir(psm_outpath)
    fprintf(1,'can not create: %s\n',psm_outpath);
    return;
end;
out_file1 = fullfile(psm_outpath,[out_filename,'.mat']);
out_file2 = fullfile(psm_outpath,[out_filename,'.plabel']);

[npep,nchg] = size(His.pep_mz);
ix = find(pep_rts(1:npep,1)>4);
if 1==isempty(ix)
    return;
end;

% terminus
nplot = length(ix);
terminus = zeros([nplot,2]);
for ino=1:nplot
    cno = ix(ino);
    p = find(isorts<=pep_rts(cno,1));
    c_ms1pos = p(end);
    c_mono_isointens = mono_isointens(:,cno);
    if pep_intens(cno,1)>0
        [nt,nb] = GetTopBottom(c_mono_isointens);%#ok
        [localmax_rt,localmax_inten,IX] = GetLocal(c_ms1pos,isorts,c_mono_isointens,nb);%#ok
        if 1==isempty(IX)
            terminus(ino,1:2) = [c_ms1pos c_ms1pos];
        else
            terminus(ino,1:2) = [IX(1) IX(end)];
        end;
    else
        terminus(ino,1:2) = [c_ms1pos c_ms1pos];
    end;
end;

% get precursors
[layoutpath,cur_raw] = fileparts(fileparts(cur_outpath));
datapath = fileparts(layoutpath);
ms2_path = fullfile(datapath,'MS2');
p0 = strfind(cur_raw,'_');
cur_rawname = cur_raw(p0(1)+1:end);
p = strfind(out_filename,'_');
c_prot = out_filename(1:p(1)-1);
sets = [0 1 2];
num_MS2 = size(MS2_index,1);
num_MS1 = size(MS1_index,1);
index = [1;MS1_index(1:num_MS1,3)];

psm.pep = repmat({''},[npep,1]);
for ino=1:npep
    psm.pep{ino,1} = [His.pep_seq,'_',His.mod_short{ino}];
end;
psm.nms2 = zeros([npep,1]);

psm.fname = {};
psm.prot = {};
psm.seq = {};
psm.mod0 = {};
psm.mod1 = {};
psm.mod2 = {};
psm.emz0 = [];
psm.emz = [];
psm.tmz = [];
psm.chg = [];
psm.rt = [];
fno = 0;
for ino=1:nplot
    if terminus(ino,1)==terminus(ino,2)
        continue;
    end;
    rt1 = isorts(terminus(ino,1));
    rt2 = isorts(terminus(ino,2));
    p = find( MS2_index(:,2)>=rt1 );
    if 1==isempty(p)
        continue;
    end;
    i1 = p(1);
    pp = find( MS2_index(:,2)<=rt2 );
    if 1==isempty(pp)
        continue;
    end;
    i2 = pp(end);

    % check MS2
    cno = ix(ino);
    flag = zeros([num_MS2,1]);
    for i=i1:i2
        cen_mz = MS2_index(i,4);
        cen_ch = MS2_index(i,5);
        for jno=1:nchg
            c_mz = His.pep_mz(cno,jno);
            c_ch = His.pep_ch(cno,jno);
            if c_ch~=cen_ch
                continue;
            end;
            mzs = c_mz + sets*unitdiff/c_ch;
            IX = find(abs(mzs-cen_mz)<ptol*cen_mz*1e-6);%#ok
            if 0==isempty(IX)
                flag(i) = jno;
                break;
            end;
        end;
    end;

    % get the info
    ms2pos = find(flag>0);
    if 0==isempty(ms2pos)
        psm.nms2(cno,1) = length(ms2pos);
        for no=1:length(ms2pos)
            i = ms2pos(no);
            cur_scan = MS2_index(i,3);
            cur_chg = MS2_index(i,5);
            fno = fno + 1;
            psm.fname{fno,1} = [cur_rawname,'.',num2str(cur_scan),'.',num2str(cur_scan),'.',num2str(cur_chg)];
            psm.prot{fno,1} = c_prot;
            psm.seq{fno,1} = His.pep_seq;
            psm.mod0{fno,1} = His.mod_short{cno};
            psm.mod1{fno,1} = His.mod_type{cno};
            psm.mod2{fno,1} = change_modtype(His.pep_seq,His.mod_type{cno});
            psm.emz0(fno,1) = MS2_index(i,4);
            psm.emz(fno,1) = MS2_index(i,4);
            c_tmz = His.pep_mz(cno,flag(i));
            psm.tmz(fno,1) = c_tmz;
            if abs(psm.emz(fno,1)-c_tmz)>ptol*c_tmz*1e-6
                p = find( MS1_index(1:num_MS1,1)<=cur_scan );
                cur_ms1pos = p(end);
                IX = index(cur_ms1pos):index(cur_ms1pos+1)-1;
                mz = MS1_peaks(IX,1);
                inten = MS1_peaks(IX,2);
                cur_ptol = ptol*c_tmz*1e-6;
                left = c_tmz - cur_ptol;
                right = c_tmz + cur_ptol;
                x = find( mz>=left & mz<=right );
                if 0==isempty(x)
                    [tmp,xx] = max(inten(x));%#ok
                    psm.emz(fno,1) = mz(x(xx));
                end;
            end;
            psm.chg(fno,1) = cur_chg;
            psm.rt(fno,1) = pep_rts(cno,1);
        end;
    end;
end;
save(out_file1,'psm');

output_plabel(out_file2,ms2_path,cur_rawname,psm);

function mod2 = change_modtype(pep,mod1)
%%

mod2 = '';

pos1 = strfind(mod1,',');
if 1==isempty(pos1);
    return;
end;

Mods = GetMods();
pos2 = [0 strfind(mod1,';')];
for jno = 1 : length(pos1)
    cpos = str2num( mod1(pos2(jno)+1:pos1(jno)-1) );%#ok
    cmod = mod1(pos1(jno)+1:pos2(jno+1)-1);
    if 0==cpos
        if 1==strcmp(cmod,'pr')%++++
          mod2 = [mod2,' 0,1'];%#ok
        end;
    else
        bflag = 0;
        no = 1;
        for i=1:length(Mods.name)
            for j=1:length(Mods.set{i})
                no = no + 1;
                % if 1==strcmp(cmod,'pr') && pep(cpos)=='K'
                if 1==strcmp(cmod,Mods.name{i}) && pep(cpos)==Mods.set{i}(j)
                    bflag = 1;
                    break;
                end;
            end;
            if 1==bflag
                break;
            end;
        end;
        mod2 = [mod2,' ',num2str(cpos),',',num2str(no)];%#ok
    end;
end;

mod2 = mod2(2:end);

function output_plabel(out_file2,ms2_path,cur_rawname,psm)
%%

fp2 = fopen(out_file2,'w');
if -1==fp2
    disp(['can not open the file: ',out_file2]);
    return;
end;

fprintf(fp2,'[FilePath]\r\n');
fprintf(fp2,'File_Path=%s\r\n',fullfile(ms2_path,[cur_rawname,'.ms2']));

Mods = GetMods();
no = 1;
fprintf(fp2,'[Modification]\r\n');
fprintf(fp2,'1=pr[PEP_N]\r\n');%++++
for i=1:length(Mods.name)
    for j=1:length(Mods.set{i})
        no = no + 1;
        % fprintf(fp2,'2=pr[K]\r\n');
        fprintf(fp2,'%d=%s[%s]\r\n',no,lower(Mods.name{i}),Mods.set{i}(j));
    end;
end;

fprintf(fp2,'[xlink]\r\n');
fprintf(fp2,'xlink=NULL\r\n');

fprintf(fp2,'[Total]\r\n');
fprintf(fp2,'total=%d\r\n',sum(psm.nms2));

% write
for ino =1 : sum(psm.nms2)
    fprintf(fp2,'[Spectrum%d]\r\n',ino);
    fname = upper(psm.fname{ino,1});
    fprintf(fp2,'name=%s\r\n',fname);
    fprintf(fp2,'pep1=0 %s 1 %s\r\n',psm.seq{ino,1},psm.mod2{ino,1});
end;

fclose(fp2);