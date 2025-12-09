function [bOK,raw_path,norganism,nsource,nsubtype] = ReadInput(para_txtfile)
%% ========================================================================
%  ReadInput
%  ------------------------------------------------------------------------
%  PURPOSE (plain English, non-proteomics audience):
%    Read a small plain-text parameter file and extract:
%      - raw_path  : folder where RAW (or general input) files live
%      - norganism : numeric code for organism (1=Human, 2=Mouse)
%      - nsource   : numeric code for data source/type
%                    (1=histone_normal, 2=histone_SILAC, 3=histone_C13, 4=histone_N15)
%      - nsubtype  : numeric subtype only relevant if nsource==4 (histone_N15)
%                    where:
%                      0: N14 light Mods
%                      1: N15 light Mods
%                      2: N14 heavy Mods
%                      3: N15 heavy Mods
%                      4: 0+1   (combine 0 and 1)
%                      5: 0+3   (combine 0 and 3)
%
%    The function scans the parameter file line-by-line to locate keys
%    'raw_path', 'norganism', 'nsource', and 'nsubtype', and returns their
%    values. It uses very simple substring searches and expects a specific
%    ordering/format (see README for details).
%
%  INPUT:
%    para_txtfile : path to a plain-text file containing one key per line,
%                   with the pattern "key=value". Example:
%                     raw_path=C:\data\run01
%                     norganism=1
%                     nsource=1
%                     nsubtype=0
%
%  OUTPUTS:
%    bOK        : 1 if the file was opened and read attempt completed;
%                 0 if the file was missing or could not be opened.
%                 % NOTE: bOK does not strictly validate that all keys were
%                 found; see "Limitations" in the README.
%    raw_path   : string (as read after '=' on the 'raw_path' line)
%    norganism  : numeric code as described above
%    nsource    : numeric code as described above
%    nsubtype   : numeric code as described above (mostly for N15 data)
%
%  IMPORTANT BEHAVIOR & ASSUMPTIONS:
%    - The code uses very basic string matching. For 'raw_path' it requires
%      the line to START with "raw_path". Leading spaces or comments will
%      make it skip lines.
%    - Keys are expected in the order: raw_path → norganism → nsource → nsubtype.
%      The scanner moves forward and does not rewind.
%    - Values are taken as the substring AFTER the first '=' on that line;
%      there is no trimming of spaces or quotes.
%    - If a key is missing, the loop may break early; defaults remain in
%      place for missing numeric values.
%
%  ========================================================================

%%
% Initialize outputs and defaults
bOK = 0;
raw_path = '';   % the datapath of raw files (string)
norganism = 1;   % default organism: 1=Human, 2=Mouse
nsource   = 1;   % default source:   1=histone_normal, 2=SILAC, 3=C13, 4=N15
nsubtype  = 0;   % for N15 only: subtype code (see header)

% -----------------------------
% Basic check: does the parameter file exist on disk?
if 0==exist( para_txtfile,'file' )
    disp(['the para file does not exist: ',para_txtfile]);
    return;
end;

% -----------------------------
% Try to open the file for reading (text mode).
fp = fopen( para_txtfile,'r' );
if -1==fp
    disp(['can not open: ',para_txtfile]);
    return;
end;

% -----------------------------
% Main read loop: proceed until end-of-file (EOF).
% Strategy:
%   1) Seek a line that starts with 'raw_path', then read its value.
%   2) Seek 'norganism' somewhere on a subsequent line; read numeric value.
%   3) Seek 'nsource'    "        "             ; read numeric value.
%   4) Seek 'nsubtype'   "        "             ; read numeric value.
% Notes:
%   - For 'raw_path', we require the match at position 1 (line must start
%     with 'raw_path'). For the others, being present anywhere on the line
%     is accepted.
%   - The code reads sequentially and does not rewind; ordering matters.
while 0==feof(fp)
    str = fgetl(fp);

    % -------------------------
    % Locate 'raw_path' at the BEGINNING of a line
    pos = strfind(str,'raw_path');
    while 0==feof(fp) && (1==isempty(pos) || 1~=pos(1))
        str = fgetl(fp);
        pos = strfind(str,'raw_path');
    end;
    if 1==feof(fp)
        break; % EOF before finding 'raw_path'
    end;
    % Extract substring after '=' and store into raw_path (no trimming)
    pos = strfind(str,'=');
    raw_path = str(pos+1:end);

    % -------------------------
    % Locate 'norganism' (can be anywhere on the line)
    pos = strfind(str,'norganism');
    while 0==feof(fp) && 1==isempty(pos)
        str = fgetl(fp);
        pos = strfind(str,'norganism');
    end;
    if 1==isempty(pos)
        break; % key missing; keep default and exit main while
    end;
    % Parse the value after '=' as a number
    pos = strfind(str,'=');
    norganism = str2double( str(pos+1:end) );

    % -------------------------
    % Locate 'nsource'
    pos = strfind(str,'nsource');
    while 0==feof(fp) && 1==isempty(pos)
        str = fgetl(fp);
        pos = strfind(str,'nsource');
    end;
    if 1==isempty(pos)
        break; % key missing; keep default and exit main while
    end;
    pos = strfind(str,'=');
    nsource = str2double( str(pos+1:end) );

    % -------------------------
    % Locate 'nsubtype'
    pos = strfind(str,'nsubtype');
    while 0==feof(fp) && 1==isempty(pos)
        str = fgetl(fp);
        pos = strfind(str,'nsubtype');
    end;
    if 1==isempty(pos)
        break; % key missing; keep default and exit main while
    end;
    pos = strfind(str,'=');
    nsubtype = str2double( str(pos+1:end) );
end;
fclose(fp);

% If we reached here, reading attempts completed; mark OK (file existed,
% opened, and we ran the scanning logic). See README for caveats about
% missing keys vs. defaults.
bOK = 1;
