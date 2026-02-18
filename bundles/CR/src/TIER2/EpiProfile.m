function EpiProfile(para_txtfile)

clc;

t1 = clock;

%% paras
% read paras
if 0==nargin || 1==isempty(para_txtfile) || 0==exist(para_txtfile,'file')
    para_txtfile = 'paras.txt';
end;
[bOK,raw_path,norganism,nsource,nsubtype] = ReadInput(para_txtfile);
if 0==bOK
    return;
end;
[def_ptol,soutput,nfigure,ndebug,raw_names] = check_otherparas(raw_path);

% check the data source
data_source = {'histone_normal','histone_SILAC','histone_C13','histone_N15','histone_13CD3'};
if ~(nsource>=1 && nsource<=length(data_source))
    fprintf(1,'please input the correct data source in [1..%d]\n',length(data_source));
    for i=1:length(data_source)
        fprintf(1,'%d: %s\n',i,data_source{i});
    end;
    return;
end;

%% raws
% raw names
if 1==isempty(raw_names)
    return;
end;

% raw2ms comentar si ya tengo ms1 y ms2
fprintf(1,'convert RAW to MS1 and MS2\n');
Raw2MS(raw_path,raw_names);

% get the MS info and ptol
ptols = repmat(def_ptol,[1,length(raw_names)]);
for i=1:length(raw_names)
    fprintf(1,'%s\n',raw_names{i});
    % get the MS1 info
    ms1_file = fullfile(raw_path,'MS1',[raw_names{i},'.MS1']);
    if 0==GetMS1ScanNo(ms1_file)
        return;
    end;
    load( fullfile(raw_path,'MS1',[raw_names{i},'_MS1scans.mat']) );
    if 1==strcmp(MS1Type,'ITMS')
        ptols(i) = 1000;
    end;

    % get the MS2 info
    ms2_file = fullfile(raw_path,'MS2',[raw_names{i},'.ms2']);
    if 0==GetMS2ScanNo(ms2_file)
        return;
    end;
end;
% ptol
if length(find(ptols==def_ptol))<length(raw_names) && length(find(ptols==1000))<length(raw_names)
    fprintf(1,'mixture of FT and IT, please separate them first!\n');
    return;
end;
ptol = ptols(1);
%{
if ndebug==1 && ptol<100
    ptol = 100;
end;
%}

diary(fullfile(raw_path,'histone_logs.txt'));

%% profiles
% get profiles
special.raw_path = raw_path;
special.nsource = nsource;
special.nsubtype = nsubtype;
special.norganism = norganism;
special.soutput = soutput;
special.nfigure = nfigure;
special.ndebug = ndebug;
if 4==nsource && (0~=nsubtype && 2~=nsubtype)
    special.nhmass = 1;
else
    special.nhmass = 0;
end;
% histone_ref
DrawISOProfile0(raw_path,raw_names,ptol,special);
% histone_normal
DrawISOProfile1(raw_path,raw_names,ptol,special);
if 2==nsource
    % histone_SILAC
    fprintf(1,'\nhistone with SILAC\n');
    DrawISOProfile2(raw_path,raw_names,ptol,special);
elseif 3==nsource
    % histone_C13
    fprintf(1,'\nhistone with C13\n');
    DrawISOProfile3(raw_path,raw_names,ptol,special);
elseif 4==nsource && (4==nsubtype || 5==nsubtype)
    % histone_N15
    fprintf(1,'\nhistone with N15\n');
    DrawISOProfile4(raw_path,raw_names,ptol,special);
elseif 5==nsource
    % histone_13CD3
    fprintf(1,'\nhistone with 13CD3\n');
    DrawISOProfile5(raw_path,raw_names,ptol,special);
end;

t2 = clock;
fprintf(['\nelapsed time: ' num2str(etime(t2,t1)) 'sec(' num2str(etime(t2,t1)/60) 'min)\n']);

diary off;