function DrawISOProfile0(raw_path,raw_names,ptol,special)
%%

% check
if ~(length(raw_names)<100 && 0==special.ndebug)
    return;
end;

layout_path = fullfile(raw_path,'histone_layouts');
if 0==exist(layout_path,'dir') && 0==mkdir(layout_path)
    fprintf(1,'can not create: %s\n',layout_path);
    return;
end;
mat_file = fullfile(layout_path,'0_ref_info.mat');
if 0~=exist(mat_file,'file')
    return;
end;

fprintf(1,'get a reference...\n');
AllUnHis = init_histone0(special);

% get top 3 raw files
ntop = 3;
range = 1:length(raw_names);
if length(raw_names)>ntop
    ipos = [1 2 3 4];% H3_01_3_8, H3_02_9_17, H3_03_18_26, H3_04_27_40
    rts = zeros(length(ipos),length(raw_names));
    intens = zeros(length(ipos),length(raw_names));
    ranks = zeros(length(ipos),length(raw_names));
    score = ones(1,length(raw_names));
    for jno=1:length(raw_names)
        cur_rawname = raw_names{jno};
        MS1_scanfile = fullfile(raw_path,'MS1',[cur_rawname,'_MS1scans.mat']);
        MS1_peakfile = fullfile(raw_path,'MS1',[cur_rawname,'_MS1peaks.mat']);
        load(MS1_scanfile);% MS1_index
        load(MS1_peakfile);% MS1_peaks
        for ino=1:length(ipos)
            [rts(ino,jno),intens(ino,jno)] = get_rts0(MS1_index,MS1_peaks,ptol,AllUnHis,ipos(ino));
        end;
    end;
    for ino=1:length(ipos)
        [tmp,xx] = sort(intens(ino,:));%#ok
        [tmp,ranks(ino,:)] = sort(xx);%#ok
        ii = find(intens(ino,:)==0);
        if 0==isempty(ii)
            ranks(ino,ii) = 1;
        end;
    end;
    for jno=1:length(raw_names)
        for ino=1:length(ipos)
            score(jno) = score(jno)*ranks(ino,jno);
        end;
    end;
    [tmp,ix] = sort(score,'descend');%#ok
    q = ix(1);
    [tmp,p] = max(intens(:,q));%#ok
    [tmp,i] = find( abs(rts(p,:)-rts(p,q))<2 );%#ok
    if length(i)<ntop
        range = i;
    else
        [tmp,x] = min( abs(rts(p,i)-median(rts(p,i))) );%#ok
        [tmp,l] = find(rts(p,i)<rts(p,i(x)));%#ok
        [tmp,r] = find(rts(p,i)>rts(p,i(x)));%#ok
        [tmp,ll] = max(intens(p,i(l)));%#ok
        [tmp,rr] = max(intens(p,i(r)));%#ok
        range = [i(l(ll)) i(x) i(r(rr))];
    end;
end;

% get rt_ref
info_rts = zeros(length(AllUnHis.pep_mz),length(range));
info_intens = zeros(length(AllUnHis.pep_mz),length(range));
for jno=1:length(range)
    cur_rawname = raw_names{range(jno)};
    MS1_scanfile = fullfile(raw_path,'MS1',[cur_rawname,'_MS1scans.mat']);
    MS1_peakfile = fullfile(raw_path,'MS1',[cur_rawname,'_MS1peaks.mat']);
    load(MS1_scanfile);% MS1_index
    load(MS1_peakfile);% MS1_peaks
    for ino=1:length(AllUnHis.pep_mz)
        [info_rts(ino,jno),info_intens(ino,jno)] = get_rts0(MS1_index,MS1_peaks,ptol,AllUnHis,ino);
    end;
    AllUnHis.selected_raws{jno,1} = cur_rawname;
end;
for ino=1:length(AllUnHis.pep_mz)
    %[tmp,ix] = max(info_intens(ino,:));%#ok
    %AllUnHis.rt_ref(ino,1) = info_rts(ino,ix);
    AllUnHis.rt_ref(ino,1) = median(info_rts(ino,:));
end;

save(mat_file,'AllUnHis');