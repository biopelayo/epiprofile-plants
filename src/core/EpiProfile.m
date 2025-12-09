function EpiProfile(para_txtfile)
% EpiProfile — Main entry point for building isotopic profiles and quantifying histone PTMs.
% ------------------------------------------------------------------------------
% INPUTS
%   para_txtfile : (char, optional) path to a text file with pipeline parameters.
%                  If missing/empty/non-existent, defaults to 'paras.txt'.
%
% SIDE EFFECTS (files/folders created/read)
%   - Reads parameters via ReadInput() and check_otherparas() from para_txtfile.
%   - Requires per-sample MS1 and MS2 files under:
%         <raw_path>/MS1/<rawname>.MS1           (custom text/binary format, project-specific)
%         <raw_path>/MS2/<rawname>.ms2           (MS2, project-specific export format)
%   - Expects <raw_path>/MS1/<rawname>_MS1scans.mat providing MS1Type (e.g., 'ITMS' vs FT).
%   - Creates log file: <raw_path>/histone_logs.txt (MATLAB diary).
%   - Creates figure/output directories used by DrawISOProfile*() helpers.
%
% NOTES
%   - This wrapper determines mass tolerance (ptol) and data source type (labeling scheme),
%     then dispatches to DrawISOProfile0..5 depending on nsource/nsubtype.
%   - Keep all raw files from the same instrument family (FT vs IT). Mixed runs are rejected.
%
% [Inference] Some exact formats (MS1/MS2/paras.txt) depend on the original EpiProfile repo.
%             Here we comment intent and control-flow without assuming undocumented details.
% ------------------------------------------------------------------------------

clc;                          % Clear MATLAB command window (purely cosmetic).

t1 = clock;                   % Start a wall-clock timer to report total elapsed time.

%% paras
% read paras
if 0==nargin || 1==isempty(para_txtfile) || 0==exist(para_txtfile,'file')
    % If no input argument, empty value, or path does not exist: default to 'paras.txt'
    para_txtfile = 'paras.txt';
end;

% ReadInput parses the parameter file and returns:
%   bOK        : success flag (1/0)
%   raw_path   : base directory containing MS1/MS2 and output folders
%   norganism  : organism selector (e.g., 1=human in original; in PLANTS, map to AT/MP/CR)
%   nsource    : data source / labeling type selector (see data_source list below)
%   nsubtype   : subtype for certain labels (e.g., special handling for 15N variants)
[bOK,raw_path,norganism,nsource,nsubtype] = ReadInput(para_txtfile);
if 0==bOK
    % If parameter parsing failed, stop early.
    return;
end;

% check_otherparas typically enriches parameters with:
%   def_ptol  : default mass tolerance (unit depends on downstream code; see below)
%   soutput   : output mode (path/name prefix or toggle, project-specific)
%   nfigure   : figure export toggle/verbosity
%   ndebug    : debug level toggle (may widen ptol in debug)
%   raw_names : cell array of sample base names to process
[def_ptol,soutput,nfigure,ndebug,raw_names] = check_otherparas(raw_path);

% check the data source
% Enumerate the supported acquisition/labeling modes.
data_source = {'histone_normal','histone_SILAC','histone_C13','histone_N15','histone_13CD3'};
if ~(nsource>=1 && nsource<=length(data_source))
    % Guardrail: if user-provided nsource is out of range, guide and exit.
    fprintf(1,'please input the correct data source in [1..%d]\n',length(data_source));
    for i=1:length(data_source)
        fprintf(1,'%d: %s\n',i,data_source{i});
    end;
    return;
end;

%% raws
% raw names
if 1==isempty(raw_names)
    % If there are no samples to process, exit gracefully.
    return;
end;

% raw2ms
fprintf(1,'convert RAW to MS1 and MS2\n');
Raw2MS(raw_path,raw_names);
% The RAW→MS1/MS2 conversion is intentionally commented out here.
% [Inference] In original pipelines, this step was done beforehand (e.g., ProteoWizard/msconvert
% + custom exporters). Uncomment and implement Raw2MS if you want in-process conversion.

% get the MS info and ptol
% Initialize per-sample mass tolerance array to the default value.
ptols = repmat(def_ptol,[1,length(raw_names)]);

for i=1:length(raw_names)
    fprintf(1,'%s\n',raw_names{i});  % Progress message: print sample name.

    % get the MS1 info
    % Build the expected MS1 file path for this sample:
    ms1_file = fullfile(raw_path,'MS1',[raw_names{i},'.MS1']);

    % Validate that the MS1 file is readable and indexable:
    % GetMS1ScanNo typically creates <rawname>_MS1scans.mat with MS1Type and index vectors.
    if 0==GetMS1ScanNo(ms1_file)
        % If MS1 parsing/indexing failed, abort.
        return;
    end;

    % Load the MS1 scans metadata produced by GetMS1ScanNo:
    % This file is expected to define a variable 'MS1Type' among other things.
    load( fullfile(raw_path,'MS1',[raw_names{i},'_MS1scans.mat']) );

    % If instrument type is Ion Trap MS (ITMS), set a coarser mass tolerance.
    % [Inference] Here ptol=1000 likely denotes a coarser window (e.g., in millimass units)
    % compared to FT instruments (Orbitrap) which use a tighter def_ptol.
    if 1==strcmp(MS1Type,'ITMS')
        ptols(i) = 1000;
    end;

    % get the MS2 info
    % Build the expected MS2 file path for this sample:
    ms2_file = fullfile(raw_path,'MS2',[raw_names{i},'.ms2']);

    % Validate MS2 file indexing/availability:
    if 0==GetMS2ScanNo(ms2_file)
        % Abort if MS2 parsing/indexing failed.
        return;
    end;
end;

% ptol
% Enforce instrument homogeneity across all samples:
% Either all runs are FT-like (all ptols==def_ptol) OR all runs are IT-like (all==1000).
if length(find(ptols==def_ptol))<length(raw_names) && length(find(ptols==1000))<length(raw_names)
    fprintf(1,'mixture of FT and IT, please separate them first!\n');
    return;
end;

% Choose the common ptol for downstream processing:
ptol = ptols(1);

%{
% Optional override while debugging:
if ndebug==1 && ptol<100
    ptol = 100;
end;
%}

% Start logging console output to a persistent text file:
diary(fullfile(raw_path,'histone_logs.txt'));

%% profiles
% get profiles
% Build the 'special' struct that carries context into DrawISOProfile* helpers.
special.raw_path  = raw_path;
special.nsource   = nsource;
special.nsubtype  = nsubtype;
special.norganism = norganism;
special.soutput   = soutput;
special.nfigure   = nfigure;
special.ndebug    = ndebug;

% Configure labeling-specific mass behavior:
% For 15N labeling (nsource==4), enable 'nhmass' unless subtype indicates exceptions.
% [Inference] 'nhmass' likely toggles heavy-isotope mass adjustments for 15N-labeled peptides.
if 4==nsource && (0~=nsubtype && 2~=nsubtype)
    special.nhmass = 1;
else
    special.nhmass = 0;
end;

% histone_ref — build reference isotope layouts (theoretical/anchor info).
DrawISOProfile0(raw_path,raw_names,ptol,special);

% histone_normal — process unlabeled (or base) histone profiles.
DrawISOProfile1(raw_path,raw_names,ptol,special);

% Dispatch to label-specific pipelines when needed:
if 2==nsource
    % histone_SILAC
    fprintf(1,'\nhistone with SILAC\n');
    DrawISOProfile2(raw_path,raw_names,ptol,special);
elseif 3==nsource
    % histone_C13
    fprintf(1,'\nhistone with C13\n');
    DrawISOProfile3(raw_path,raw_names,ptol,special);
elseif 4==nsource && (4==nsubtype || 5==nsubtype)
    % histone_N15 — only for certain subtypes (4 or 5).
    fprintf(1,'\nhistone with N15\n');
    DrawISOProfile4(raw_path,raw_names,ptol,special);
elseif 5==nsource
    % histone_13CD3 (e.g., propionyl-d3 derivatization or similar heavy tags).
    fprintf(1,'\nhistone with 13CD3\n');
    DrawISOProfile5(raw_path,raw_names,ptol,special);
end;

% Report elapsed time (seconds and minutes):
t2 = clock;
fprintf(['\nelapsed time: ' num2str(etime(t2,t1)) 'sec(' num2str(etime(t2,t1)/60) 'min)\n']);

% Stop logging:
diary off;
end
