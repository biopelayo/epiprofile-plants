function [def_ptol,soutput,nfigure,ndebug,raw_names] = check_otherparas(raw_path)
%% check_otherparas — run-time defaults and the list of runs (raw_names)
%
%  TIER: T2 (modified from upstream EpiProfile 2.0 basic)
%
%  raw_names is the list of run basenames. Upstream derives it from
%  <raw_path>/*.raw only. PLANTS keeps that and, when no .raw is present,
%  falls back to the basenames of <raw_path>/MS1/*.MS1 (or *.ms1), so data
%  converted from .wiff/.d/mzML no longer needs 0-byte .raw placeholders.
%  Both lists are sorted, so run numbering (01_, 02_, ...) is deterministic
%  across operating systems (Fix A-CRIT-004 / C-CRIT-04, 2026-05-24; dir()
%  order is not guaranteed on Windows).

% **default can be changed**
def_ptol = 10;% FT: 10 ppm, or others < 100 ppm
soutput = '31';% 1st bit: 1, H3H4 basic ; 2, H3H4 basic + H3S10ph; 3, H3H4 all    2nd bit: 1, H1H2AB; 0, no H1H2AB
nfigure = 0;% 1: output figures (needs Statistics + Bioinformatics Toolboxes and a display;
            %    print -dpdf failed under matlab -batch, F4 2026-06-16), 0: no figures
ndebug = 0;% 0: normal and ref, 1: debug, 2: normal but no ref

% raw_path\raw_names
raws = [dir(fullfile(raw_path,'*.raw')); dir(fullfile(raw_path,'*.RAW'))];
if 1==isempty(raws)
    ms1s = [dir(fullfile(raw_path,'MS1','*.MS1')); dir(fullfile(raw_path,'MS1','*.ms1'))];
    if 1==isempty(ms1s)
        raw_names = {};
        fprintf(1,'no raws: neither *.raw in %s nor *.MS1 in its MS1/ subfolder\n',raw_path);
        return;
    end;
    names = cell(length(ms1s),1);
    for i=1:length(ms1s)
        [~,names{i}] = fileparts(ms1s(i).name);
    end;
    raw_names = unique(names);% unique() sorts; a case-insensitive filesystem lists the same file twice
    raw_names = raw_names(:);
    fprintf(1,'no *.raw found, run list taken from MS1/ (%d runs)\n',length(raw_names));
    return;
end;

names = cell(length(raws),1);
for i=1:length(raws)
    names{i} = raws(i).name(1:end-4);
end;
raw_names = unique(names);% sorted and de-duplicated (.raw/.RAW on case-insensitive filesystems)
raw_names = raw_names(:);
