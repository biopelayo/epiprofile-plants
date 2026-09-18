function EpiProfile(para_txtfile)
%% EpiProfile — entry point of an EpiProfile_PLANTS species bundle
%
%  TIER: T2 (modified from upstream EpiProfile 2.0 basic)
%
%  USAGE
%    EpiProfile                 % reads ./paras.txt (current folder)
%    EpiProfile('C:\data\paras.txt')
%
%  paras.txt (see paras.example.txt at the repository root):
%    [EpiProfile]
%    raw_path=<folder with <run>.raw, MS1/ and MS2/>
%    norganism=1
%    nsource=1        % histone_normal is the only mode ported to PLANTS
%    nsubtype=0
%
%  Expected data layout under raw_path (MS1/MS2 already extracted):
%    <raw_path>/<run>.raw          names only; 0-byte placeholders are fine
%    <raw_path>/MS1/<run>.MS1      (.ms1 also accepted)
%    <raw_path>/MS2/<run>.ms2      (.MS2 also accepted)
%  If MS1/MS2 are missing, RawToMS1.exe and xtract.exe must sit in the
%  current folder (pwd) and <run>.RAW must be a real Thermo file (Raw2MS.m).
%
%  Changes vs upstream:
%   2026-05-24 (Tribunal de Codigo, capas A/G y C-CRIT-05): no clc (it hung
%     MATLAB in -batch mode); startup timestamp so an empty log means "hung
%     before parsing"; abort if DrawISOProfile0 did not write 0_ref_info.mat
%     (otherwise DrawISOProfile1 continues on a failed calibration).
%   2026-08-18: only nsource==1 accepted (the SILAC/C13/N15/13CD3 runners were
%     not ported), Raw2MS only invoked when some MS1/MS2 file is missing,
%     MS1/MS2 extensions matched case-insensitively, clearer paras messages.

fprintf(1,'EpiProfile started: %s\n', datestr(now));
t1 = clock;

%% paras
% read paras
if 0==nargin || 1==isempty(para_txtfile)
    para_txtfile = 'paras.txt';
end;
if 0==exist(para_txtfile,'file')
    fprintf(1,'paras file not found: %s (pwd: %s)\n',para_txtfile,pwd);
    fprintf(1,'Pass its path, EpiProfile(''C:\\data\\paras.txt''), or run from a folder that contains paras.txt.\n');
    fprintf(1,'A template lives at the repository root: paras.example.txt\n');
    return;
end;
[bOK,raw_path,norganism,nsource,nsubtype] = ReadInput(para_txtfile);
if 0==bOK
    return;
end;
if 0==exist(raw_path,'dir')
    fprintf(1,'raw_path does not exist: %s (check paras.txt)\n',raw_path);
    return;
end;
[def_ptol,soutput,nfigure,ndebug,raw_names] = check_otherparas(raw_path);

% check the data source
% Only histone_normal (nsource=1) is ported to the PLANTS bundles. The
% upstream runners for the other modes (DrawISOProfile3/4/5, Extract_SILAC)
% are not in the bundle, so fail here with a clear message instead of an
% "Undefined function" error after the (long) MS1/MS2 parsing.
data_source = {'histone_normal','histone_SILAC','histone_C13','histone_N15','histone_13CD3'};
if 1~=nsource
    fprintf(1,'nsource=%d is not supported by EpiProfile_PLANTS; only 1 (%s) is ported.\n',nsource,data_source{1});
    fprintf(1,'Upstream modes not ported:\n');
    for i=2:length(data_source)
        fprintf(1,'%d: %s\n',i,data_source{i});
    end;
    return;
end;

%% raws
% raw names
if 1==isempty(raw_names)
    return;
end;

% raw2ms: only when some MS1/MS2 file is missing (needs RawToMS1.exe and
% xtract.exe in pwd, see Raw2MS.m). With pre-extracted MS1/MS2 nothing is run.
need_convert = false;
for i=1:length(raw_names)
    if isempty(find_ms_file(fullfile(raw_path,'MS1'),raw_names{i},{'.MS1','.ms1'})) || ...
       isempty(find_ms_file(fullfile(raw_path,'MS2'),raw_names{i},{'.ms2','.MS2'}))
        need_convert = true;
        break;
    end;
end;
if need_convert
    fprintf(1,'convert RAW to MS1 and MS2\n');
    Raw2MS(raw_path,raw_names);
else
    fprintf(1,'MS1/MS2 files already present, skipping RAW conversion\n');
end;

% get the MS info and ptol
ptols = repmat(def_ptol,[1,length(raw_names)]);
for i=1:length(raw_names)
    fprintf(1,'%s\n',raw_names{i});
    % get the MS1 info
    ms1_file = find_ms_file(fullfile(raw_path,'MS1'),raw_names{i},{'.MS1','.ms1'});
    if isempty(ms1_file)
        ms1_file = fullfile(raw_path,'MS1',[raw_names{i},'.MS1']);% for the error message
    end;
    if 0==GetMS1ScanNo(ms1_file)
        return;
    end;
    load( fullfile(raw_path,'MS1',[raw_names{i},'_MS1scans.mat']) );
    if 1==strcmp(MS1Type,'ITMS')
        ptols(i) = 1000;
    end;

    % get the MS2 info
    ms2_file = find_ms_file(fullfile(raw_path,'MS2'),raw_names{i},{'.ms2','.MS2'});
    if isempty(ms2_file)
        ms2_file = fullfile(raw_path,'MS2',[raw_names{i},'.ms2']);
    end;
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
% Fix C-CRIT-05 (2026-05-24): make sure the calibration wrote its reference.
% Without this, a silent failure in DrawISOProfile0 let the pipeline continue
% and produce corrupt outputs (seen on PXD014739 and PXD046034).
ref_mat_file = fullfile(raw_path,'histone_layouts','0_ref_info.mat');
if 2~=exist(ref_mat_file,'file')
    fprintf(1,'ERROR: calibration failed (0_ref_info.mat was not created). Aborting.\n');
    diary off;
    return;
end;
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

%--------------------------------------------------------------------------
function f = find_ms_file(folder,basename,exts)
% first existing <folder>/<basename><ext> for the given extensions, '' if none.
% Keeps the upstream on-disk convention (.MS1 / .ms2) but also accepts the
% other case, which matters on case-sensitive filesystems.
f = '';
for k=1:length(exts)
    cand = fullfile(folder,[basename,exts{k}]);
    if 0~=exist(cand,'file')
        f = cand;
        return;
    end;
end;
