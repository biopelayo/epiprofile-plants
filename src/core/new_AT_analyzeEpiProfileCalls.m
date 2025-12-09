% analyzeEpiProfileCalls.m
% This utility scans a MATLAB folder, locates `EpiProfile.m`, extracts the names
% of (sub)functions declared inside that file, runs MATLAB's profiler while
% executing `EpiProfile`, and then prints a color-coded list in the console:
%   - Green: subfunction that was actually called during the profiled run
%   - Red  : subfunction that was never called
%
% Usage:
%   analyzeEpiProfileCalls('path/to/folder/with/EpiProfile.m');

function analyzeEpiProfileCalls(dirPath)
    % ANALYZEEPIPROFILECALLS — Static subfunction discovery + dynamic call check
    %
    % INPUT
    %   dirPath (char/string): directory containing EpiProfile.m
    %
    % SIDE EFFECTS
    %   - Adds dirPath to the MATLAB path (temporarily for this session).
    %   - Starts/stops MATLAB profiler and queries profile info.
    %   - Prints a colored list of subfunctions to the Command Window
    %     using ANSI escape sequences (see README notes).
    %
    % LIMITATIONS (by design, logic intentionally unchanged)
    %   - Subfunction discovery is purely textual (regex), based solely on
    %     the contents of EpiProfile.m. It does NOT parse/inspect other files.
    %   - The name comparison between "declared subfunctions" and "executed
    %     entries from the profiler" is string-based and exact. If MATLAB
    %     reports names with scope/package/class prefixes, they may not match
    %     plain subfunction names and therefore appear as "never called"
    %     even if they executed.
    %   - The code assumes `EpiProfile` can run with default/no arguments.
    %     If it requires inputs or setup, you will need to provide those
    %     within EpiProfile itself (logic not modified here).
    %
    % NOTE ON COLORS
    %   - The printout uses ANSI escape codes (\x1b[32m, \x1b[31m). Some MATLAB
    %     Command Window environments do not render ANSI colors; in those cases,
    %     you will just see the raw sequences. Running MATLAB in a system terminal
    %     (e.g., -nodesktop) may render colors depending on OS/terminal support.

    % --- 0) Basic input validation ------------------------------------------------
    if nargin < 1 || ~isfolder(dirPath)
        error('Debe proporcionar una carpeta válida donde están los scripts MATLAB.');
    end

    % --- 1) Locate the main script file: EpiProfile.m -----------------------------
    mainScript = fullfile(dirPath, 'EpiProfile.m');
    if ~isfile(mainScript)
        error('No se encontró EpiProfile.m en %s', dirPath);
    end

    % --- 2) Read the file and extract subfunction names via regex -----------------
    % We read the entire file as text, then iterate line-by-line looking for
    % "function <name>(" or "function [out,...] = <name>(" patterns.
    text = fileread(mainScript);
    lines = splitlines(text);
    subfuncs = {};
    % Start at line 2 to intentionally skip the primary function declaration
    % if it appears on the very first line. This preserves the original logic.
    for i = 2:numel(lines)
        % REGEX EXPLANATION:
        % ^\s*function                 : line starting with (optional) whitespace + 'function'
        % (?:\[[^\]]*\]\s*=)?          : optional output pattern in brackets `[ ... ] =`
        % \s*(\w+)\s*\(                : capture the function name (word chars) before '('
        tok = regexp(lines{i}, ...
            '^\s*function\s+(?:\[[^\]]*\]\s*=)?\s*(\w+)\s*\(', ...
            'tokens','once');
        if ~isempty(tok)
            subfuncs{end+1} = tok{1}; %#ok<AGROW>
        end
    end
    % Remove potential duplicates (e.g., multiple declarations with same name)
    subfuncs = unique(subfuncs);
    if isempty(subfuncs)
        fprintf('No se encontraron subfunciones en EpiProfile.m\n');
        return;
    end

    % --- 3) Profile the execution of EpiProfile ----------------------------------
    % We reset the profiler, add the directory to path, attempt to run EpiProfile,
    % and then stop the profiler to query call information.
    profile off; profile clear; profile on -timer real;
    addpath(dirPath);
    try
        % NOTE: If EpiProfile requires inputs, you would adapt the call here.
        % Logic intentionally left unchanged as requested.
        EpiProfile;
    catch ME
        % We do not abort; we still want to list any functions that were called
        % before the error occurred. Hence, warn and proceed.
        warning('Error al ejecutar EpiProfile: %s', ME.message);
    end
    profile off;

    % --- 4) Inspect profiler info -------------------------------------------------
    % 'profile(''info'')' returns a struct with a FunctionTable containing entries
    % for functions that were actually observed during the profiling session.
    pInfo = profile('info');
    executed = unique({pInfo.FunctionTable.FunctionName});

    % --- 5) Pretty-print results (green = executed, red = never called) ----------
    % IMPORTANT: We compare plain names from `subfuncs` with the profiler's
    % FunctionName strings. If profiler reports names with scope qualifiers
    % (e.g., packages/classes), exact matching may fail and the line will
    % print in red even if something similar was executed.
    fprintf('\nResultados de ejecución de EpiProfile:\n');
    for k = 1:numel(subfuncs)
        name = subfuncs{k};
        if any(strcmp(name, executed))
            % GREEN (ANSI 32)
            fprintf('\x1b[32m%s\x1b[0m\n', name);
        else
            % RED (ANSI 31)
            fprintf('\x1b[31m%s\x1b[0m\n', name);
        end
    end
end
