function check_layout()
% the cursor goes to the last line, and press F12 to set a breakpoint
% press F5 to start, and check layouts for RAW
% if saving changes, press F5; otherwise, press red button above 'Quit Debugging'
%
% CHECK_LAYOUT  Interactive helper to load histone layout metadata, inspect/edit it, and optionally save.
%
% Synopsis:
%   check_layout()
%
% Purpose:
%   This routine implements a minimal, *interactive* workflow to review and (optionally) persist
%   changes to reference information stored in a MAT-file named:
%       <raw_path>/histone_layouts/0_ref_info.mat
%   The intended usage is to run the function under the MATLAB debugger with a breakpoint on
%   the **last line** (the SAVE call). While paused *before* saving, you can inspect and edit
%   variables (especially `AllUnHis`) in the function workspace. Then either:
%     - Press **F5** (Continue) to execute the final `save(...)` and persist changes to disk, or
%     - Click **Quit Debugging** (red square) to stop execution *before* saving, discarding edits.
%
% Inputs / external resources:
%   - A text configuration file 'paras.txt' is read via `ReadInput('paras.txt')`, which is expected
%     to return a success flag and a `raw_path` directory. [Inferencia] El contenido exacto de
%     'paras.txt' y el contrato de `ReadInput` no están definidos aquí; típicamente incluirá rutas.
%
% Outputs:
%   - None returned. Side-effect only: if you choose to continue past the breakpoint, the function
%     overwrites (or creates) the MAT-file and writes variable `AllUnHis` into it.
%
% Side effects:
%   - Clears the Command Window (`clc`).
%   - Loads variables from the MAT-file into the *function workspace* (not base workspace).
%   - Optionally saves variable `AllUnHis` back to the same MAT-file.
%
% Requirements / dependencies:
%   - Function `ReadInput` must exist on path and accept 'paras.txt', returning:
%         [bOK, raw_path]
%     where `bOK` is nonzero on success and `raw_path` is a valid folder.
%   - The file `<raw_path>/histone_layouts/0_ref_info.mat` should exist and contain
%     `AllUnHis` among other metadata (this code will error on `save` if `AllUnHis` is
%     undefined in the workspace). [Inferencia]
%
% Debugging workflow (recommended):
%   1) Place the cursor on the **last line** (`save(...)`) and press **F12** to set a breakpoint.
%   2) Press **F5** to run:
%        - the function will `clc`, read `paras.txt` to obtain `raw_path`, check the MAT-file,
%          and `load` it into the function workspace.
%        - execution pauses at the last line (due to your breakpoint).
%   3) While paused, inspect/edit variables (e.g., `AllUnHis`) in the workspace.
%   4) To **save** your edits: press **F5** to continue and execute `save(...)`.
%      To **discard** edits: click **Quit Debugging** (red square) so `save(...)` is not executed.
%
% Notes / caveats:
%   - If `ReadInput` fails or the MAT-file is missing, the function returns early (no error).
%   - `load(mat_file)` without output injects all variables in the MAT-file into the function's
%     local workspace. They are not persisted unless you continue past the breakpoint to `save`.
%   - `save(mat_file,'AllUnHis')` will throw if `AllUnHis` is not present in the workspace.
%   - `exist(mat_file,'file')` returns 2 for files; the code checks equality to zero to early-return.
%

clc;
% Clear Command Window to provide a clean view for the interactive session.

[bOK,raw_path] = ReadInput('paras.txt');% read paras
% Read configuration parameters from 'paras.txt'.
% Expected: bOK (1 on success, 0 on failure), and raw_path (directory path).
% If bOK==0, we abort silently (no error), preserving a non-intrusive interactive workflow.

if 0==bOK
    return;
end;

mat_file = fullfile(raw_path,'histone_layouts','0_ref_info.mat');
% Build the absolute path to the reference MAT-file inside <raw_path>/histone_layouts/.

if 0==exist(mat_file,'file')
    return;
end;
% If the reference file does not exist, abort quietly (no error).

load(mat_file);% load ref
% Load all variables from the MAT-file into the function workspace.
% Typical content (by convention) includes 'AllUnHis' (histone layouts). [Inferencia]

fprintf(1,'checking layouts...\n');
% Informational message to the Command Window / fd=1 (stdout).

save(mat_file,'AllUnHis');% save changes
% Persist ONLY the variable 'AllUnHis' back into the same MAT-file.
% IMPORTANT: Place a breakpoint here (F12) before running. While paused, inspect/edit
% 'AllUnHis'. Then press F5 to execute this line and write changes, or Quit Debugging to skip it.
