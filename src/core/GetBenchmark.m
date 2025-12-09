function GetBenchmark(cur_outpath)
%% GetBenchmark — Export a tab-delimited identification summary from PSM .mat files
%
% PURPOSE (plain English, non-expert friendly)
%   This utility scans the subfolder "<cur_outpath>/psm" for .mat files that
%   contain a variable named 'psm' (Peptide-Spectrum Match table). For every
%   entry in each 'psm', it computes the precursor mass error in parts-per-million
%   (ppm) and writes a unified, tab-delimited report called:
%
%       <cur_outpath>/psm/identification_list.xls
%
%   NOTE: Despite the ".xls" extension, the file is a plain *tab-separated text*
%   file produced with fprintf (not a native Excel binary). Excel/LibreOffice can
%   open it directly, but it is technically a TSV.
%
% EXPECTED CONTENTS OF 'psm' (loaded from each .mat):
%   psm.fname : cell array of char, filename per PSM (e.g., raw/mzML source)
%   psm.emz   : vector double, *measured* precursor m/z (experimental)
%   psm.tmz   : vector double, *theoretical* precursor m/z (from sequence/mods)
%   psm.chg   : vector int,    precursor charge state
%   psm.seq   : cell array of char, peptide sequence (1-letter, uppercase)
%   psm.mod0  : cell array of char, modification string #1 (project-specific)
%   psm.mod1  : cell array of char, modification string #2 (project-specific)
%   psm.prot  : cell array of char, protein/histone type annotation
%   psm.rt    : vector double, retention time (minutes)
%
% OUTPUT
%   A single TSV-like file with header and one row per PSM across all .mat inputs.
%
% COLUMN LAYOUT OF THE OUTPUT FILE
%   1) filename
%   2) measured m/z
%   3) calculated m/z
%   4) charge
%   5) ppm
%   6) sequence
%   7) modification1
%   8) modification2
%   9) histone type
%  10) retention time(min)
%
% PPM FORMULA
%   ppm = 1e6 * (emz - tmz) / tmz
%   (positive = measured is larger than theoretical)
%
% CAVEATS
%   - The code assumes each .mat contains a variable named 'psm' with the fields
%     listed above. If not present, 'load' will still succeed but the subsequent
%     field access will error. Validate your inputs if needed.
%   - If any text field contains TABs or newlines, the TSV alignment can break.
%     Keep text fields simple (no tabs/newlines) or sanitize beforehand.
%   - The output is opened with default text mode. If you work on non-UTF-8
%     environments or need consistent line endings, adapt fopen options accordingly.
%
% USAGE
%   GetBenchmark('/path/to/output_root');
%   % Expects '/path/to/output_root/psm/*.mat' to exist.
%

    % Derive the PSM folder inside the provided output path
    psm_outpath = fullfile(cur_outpath,'psm');

    % Compute target report file path and attempt to open it for writing
    benchmarkfile = fullfile(psm_outpath,'identification_list.xls');
    fp = fopen(benchmarkfile,'w');   % plain text, write mode
    if -1 == fp
        % If opening fails (e.g., missing permissions), inform and exit gracefully
        fprintf('can not open:%s\n', benchmarkfile);
        return;
    end

    % Write the header line (TAB-delimited). This defines the column order.
    fprintf(fp, ['filename\tmeasured m/z\tcalculated m/z\tcharge\tppm\t', ...
                 'sequence\tmodification1\tmodification2\thistone type\t', ...
                 'retention time(min)\n']);

    % List all .mat files under <cur_outpath>/psm
    matfiles = dir(fullfile(psm_outpath,'*.mat'));

    % Iterate over each .mat file
    for ino = 1:length(matfiles)
        % Compose absolute path and load variables into the current workspace
        cur_file = fullfile(psm_outpath, matfiles(ino).name);
        load(cur_file);  %#ok<LOAD>  % expects a variable named 'psm'

        % Number of PSM entries (length of fname defines the row count)
        nlen = length(psm.fname);

        % Loop over all rows in the 'psm' table
        for jno = 1:nlen
            % Compute mass accuracy in ppm using measured (emz) vs theoretical (tmz)
            ppm = 1e6 * (psm.emz(jno) - psm.tmz(jno)) / psm.tmz(jno);

            % Emit one tab-separated row per PSM with controlled numeric formats:
            %   m/z to 6 decimals, ppm to 6 decimals, RT to 4 decimals, charge as integer
            fprintf(fp, '%s\t%.6f\t%.6f\t%d\t%.6f\t%s\t%s\t%s\t%s\t%.4f\n', ...
                psm.fname{jno}, psm.emz(jno), psm.tmz(jno), psm.chg(jno), ppm, ...
                psm.seq{jno}, psm.mod0{jno}, psm.mod1{jno}, psm.prot{jno}, psm.rt(jno));
        end
    end

    % Always close the file handle to flush buffers and release OS resources
    fclose(fp);
end
