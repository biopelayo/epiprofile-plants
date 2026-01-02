function [def_ptol,soutput,nfigure,ndebug,raw_names] = check_otherparas(raw_path)
%%
% CHECK_OTHERPARAS  Prepare default processing parameters and enumerate RAW files in a folder.
%
% Syntax:
%   [def_ptol, soutput, nfigure, ndebug, raw_names] = check_otherparas(raw_path)
%
% Purpose:
%   This helper initializes a set of *default* analysis/display parameters and collects the list
%   of RAW file basenames present under the provided directory. It is typically used early in a
%   workflow to standardize tolerances/flags and to know which runs (RAWs) are available.
%
% Input:
%   - raw_path : char/string. Absolute or relative path to the folder containing Thermo *.RAW files.
%                Note: the code looks for the **uppercase** extension '*.RAW'. On case-sensitive
%                filesystems, lowercase '.raw' will not match.
%
% Outputs:
%   - def_ptol : double. Default mass tolerance in parts-per-million (ppm). The in-line note says
%                “FT: 10 ppm, or others < 100 ppm”, suggesting this is meant for high-resolution data.
%   - soutput  : char (string). Two-character code controlling which outputs are produced:
%                  * 1st digit (H3/H4 scope): '1' = H3H4 basic; '2' = H3H4 basic + H3S10ph; '3' = H3H4 all
%                  * 2nd digit (H1/H2AB scope): '1' = include H1/H2AB; '0' = exclude H1/H2AB
%                Example: '30' => H3/H4 = all; H1/H2AB = no.
%   - nfigure  : double (0/1). Flag to control figure generation:
%                  1 = produce figures; 0 = no figures.
%   - ndebug   : double (mode). Debugging/reference behavior:
%                  0 = normal + use reference; 1 = debug; 2 = normal but **no** reference.
%                (Exact downstream effects depend on consumer functions.)
%   - raw_names: cellstr (N×1). Basenames (without extension) of all '*.RAW' files found in raw_path.
%                If none are found, returns {} and prints 'no raws'.
%
% Behavior:
%   1) Sets the defaults (def_ptol, soutput, nfigure, ndebug).
%   2) Lists files matching '<raw_path>/*.RAW'.
%   3) If none found: returns with raw_names = {}, printing a notice.
%   4) Otherwise: returns the list of RAW basenames (without '.RAW') in raw_names.
%
% Dependencies:
%   - Built-ins: dir, fullfile, repmat, fprintf.
%   - No project-specific functions are called here.
%
% Notes / Caveats:
%   - Case sensitivity: on Linux/macOS with case-sensitive FS, files named '*.raw' won’t match.
%   - Path validity: if raw_path is invalid, dir(...) returns empty, leading to 'no raws'.
%   - soutput is a **string**, not a numeric mask; downstream code must parse it accordingly.
%

% **default can be changed**
def_ptol = 10;% FT: 10 ppm, or others < 100 ppm
% Default monoisotopic mass tolerance in ppm (typically for FT/Orbitrap-level HR-MS).

soutput = '10';% 1st bit: 1, H3H4 basic ; 2, H3H4 basic + H3S10ph; 3, H3H4 all    2nd bit: 1, H1H2AB; 0, no H1H2AB
% Two-character output-control code as string:
%   - First character controls H3/H4 scope (1=basic, 2=basic+H3S10ph, 3=all).
%   - Second character toggles H1/H2AB inclusion (1=include, 0=exclude).

nfigure = 1;% 1: output figures, 0: no figures
% Whether to generate figures downstream (consumer-dependent).

ndebug = 0;% 0: normal and ref, 1: debug, 2: normal but no ref
% Debug/Reference mode selector. Actual meaning is interpreted by later stages.

% raw_path\raw_names
raws = dir(fullfile(raw_path,'*.RAW'));
% Enumerate all files ending with '.RAW' under raw_path.
% Note: On case-sensitive systems, only uppercase '.RAW' will match.

if 1==isempty(raws)
    raw_names = {};
    fprintf(1,'no raws\n');
    return;
end;
% If no RAW files are found, return early with an empty list and a console message.

raw_names = repmat({''},[length(raws),1]);
% Preallocate a cell array (N×1) to hold basenames for each RAW file found.

for i=1:length(raws)
    raw_names{i,1} = raws(i).name(1:end-4);
end;
% Strip the last four characters ('.RAW') to keep only the basename for each entry.
