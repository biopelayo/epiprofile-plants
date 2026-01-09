function [def_ptol,soutput,nfigure,ndebug,raw_names] = check_otherparas(raw_path)
%%

% **default can be changed**
def_ptol = 10;% FT: 10 ppm, or others < 100 ppm
soutput = '11';% 1st bit: 1, H3H4 basic ; 2, H3H4 basic + H3S10ph; 3, H3H4 all    2nd bit: 1, H1H2AB; 0, no H1H2AB
nfigure = 1;% 1: output figures, 0: no figures
ndebug = 0;% 0: normal and ref, 1: debug, 2: normal but no ref

% raw_path\raw_names
raws = dir(fullfile(raw_path,'*.RAW'));
if 1==isempty(raws)
    raw_names = {};
    fprintf(1,'no raws\n');
    return;
end;

raw_names = repmat({''},[length(raws),1]);

for i=1:length(raws)
    raw_names{i,1} = raws(i).name(1:end-4);
end;