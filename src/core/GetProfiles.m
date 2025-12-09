function [isorts,c_ref_isointens] = GetProfiles(MS1_index,MS1_peaks,c_ref_isomzs,cur_ch,ptol,nC13,PosRange)
%% GetProfiles — Extract MS1 chromatographic traces at reference isotopic m/z
%
% PLAIN-ENGLISH PURPOSE (for non-experts)
%   For a target precursor (known m/z and charge), we want the MS1
%   **extracted ion chromatograms (XICs)** at the expected isotopic
%   positions (e.g., M−1, M, M+1, M+2...). This function:
%     1) walks scan-by-scan over a given subset of scans (PosRange),
%     2) around each reference m/z window, finds candidate peaks, and
%     3) stores per-scan intensities for each isotopic channel.
%
%   Two operation modes:
%     • nC13 == 1  -> SIMPLE mode: take the maximum intensity found within
%                     each m/z window independently (fast, permissive).
%     • nC13 ~= 1  -> PATTERN-AWARE mode: enforce basic isotopic consistency
%                     using a precomputed table (IPV) of expected M+1/M
%                     ratios for a peptide of similar mass. This helps
%                     disambiguate overlapping signals.
%
% INPUTS
%   MS1_index       [nScans × 4]  -> columns: [scan_no, rt, ptr, baseline]
%   MS1_peaks       [nPeaks × 2]  -> concatenated [m/z, intensity] across scans
%   c_ref_isomzs    [1 × nmz]     -> reference isotopic m/z values
%                                    (e.g., [M−1, M, M+1, M+2], same order downstream)
%   cur_ch          scalar        -> precursor charge state (z)
%   ptol            scalar (ppm)  -> m/z tolerance (ppm); a hard cap of 0.3 Th is applied
%   nC13            flag          -> mode selector:
%                                     1 = SIMPLE (no isotopic ratio checks)
%                                     else = PATTERN-AWARE (use IPV table)
%   PosRange        vector of int -> MS1 scan indices to process (e.g., rt_i1:rt_i2)
%
% OUTPUTS
%   isorts          [nScans × 1]  -> retention time (RT) per scan (minutes)
%   c_ref_isointens [nScans × nmz]-> per-scan intensity for each isotopic channel
%                                    (0 if nothing acceptable was found at that scan)
%
% HOW TO READ PEAKS FOR A SCAN
%   The peaks for scan i are stored in rows:
%       IX = [ MS1_index(i,3) : MS1_index(i+1,3)-1 ]   (or up to end for last scan)
%       mz   = MS1_peaks(IX,1);
%       inten= MS1_peaks(IX,2);
%
% NOTES
%   • The isotopic reference array (c_ref_isomzs) must match the channel
%     meaning used later (e.g., j=2 is the monoisotopic M, j=3 is M+1, etc.).
%   • PATTERN-AWARE mode uses an IPV table whose row for mass ~M0 holds
%     relative intensities [I(M), I(M+1), I(M+2), I(M+3), I(M+4)] normalized
%     to I(M)=1. We only use the M+1/M ratio (t1) here.
%   • Proton mass used for neutral mass estimate: 1.007276 Da.

    % ---- 0) Setup and preallocation ----
    nmz     = length(c_ref_isomzs);              % number of isotopic channels
    num_MS1 = size(MS1_index,1);

    % 'index' gives the start pointer in MS1_peaks for each scan
    % index(i) .. index(i+1)-1 are the rows for MS1 scan i
    index = [1; MS1_index(1:num_MS1,3)];

    % RT per scan (copied out for convenience)
    isorts = MS1_index(1:num_MS1,2);

    % Output matrix: one column per isotopic m/z, one row per scan
    c_ref_isointens = zeros([num_MS1, nmz]);

    % ---- 1) SIMPLE mode: independent windowed maxima (fast & permissive) ----
    if 1==nC13
        for ino = PosRange
            % Slice the peaks of MS1 scan 'ino'
            IX    = index(ino):index(ino+1)-1;
            mz    = MS1_peaks(IX,1);
            inten = MS1_peaks(IX,2);

            % For each isotopic window, take the maximum intensity inside tolerance
            tmp_intens = zeros([1, nmz]);
            for jno = 1:nmz
                c_mz   = c_ref_isomzs(jno);
                % Convert ppm to absolute Th, cap at 0.3 Th to avoid huge windows
                c_ptol = min([ptol*c_mz*1e-6, 0.3]);
                left   = c_mz - c_ptol;
                right  = c_mz + c_ptol;

                pos = find( mz>=left & mz<=right );
                if ~isempty(pos)
                    % peak picking here is just "take the local max in the window"
                    tmp_intens(jno) = max(inten(pos));
                end
            end
            c_ref_isointens(ino,1:nmz) = tmp_intens;
        end

    % ---- 2) PATTERN-AWARE mode: use isotope ratios to select consistent peaks ----
    else
        thr = 5;                 % allow deviations up to ×5 from the expected ratio (log-space check)
        IPV = GetIPV();          % precomputed isotope pattern vectors (I(M)=1, then M+1..M+4 fractions)

        % Estimate nominal peptide (neutral) mass to select IPV row:
        %   neutral_mass ≈ M (monoisotopic m/z at j=2) * z − z*H+
        % We use floor(M0) to index IPV rows (as in the original code).
        M    = c_ref_isomzs(2)*cur_ch - cur_ch*1.007276;  % neutral mass estimate
        tIPV = IPV(floor(M), 1:5);                        % [I(M), I(M+1), ...], normalized to I(M)=1
        t1   = tIPV(2);                                   % expected ratio I(M+1)/I(M)

        for ino = PosRange
            % Slice peaks for this scan
            IX    = index(ino):index(ino+1)-1;
            mz    = MS1_peaks(IX,1);
            inten = MS1_peaks(IX,2);

            % 2.1 Flag positions that fall inside each isotopic window
            %     flags(k)=j means peak k is inside window of channel j
            flags = zeros([1, length(mz)]);
            for jno = 1:nmz
                c_mz   = c_ref_isomzs(jno);
                c_ptol = min([ptol*c_mz*1e-6, 0.3]);
                left   = c_mz - c_ptol;
                right  = c_mz + c_ptol;
                pos = find( mz>=left & mz<=right );
                if ~isempty(pos)
                    flags(pos) = jno;
                end
            end

            % 2.2 Gather candidate indices for M−1, M, M+1, M+2 (assuming nmz>=4)
            %     NOTE: channel order must match 'c_ref_isomzs' passed in
            tmp_intens = zeros([1, nmz]);
            x1 = find(flags==1);   % candidates in M−1 window
            x2 = find(flags==2);   % candidates in M   window (monoisotopic)
            x3 = find(flags==3);   % candidates in M+1 window
            x4 = find(flags==4);   % candidates in M+2 window

            % 2.3 Choose the (M, M+1) pair whose intensity ratio best matches 't1'
            if ~isempty(x2) && ~isempty(x3)
                % Build all pairs (x2, x3) and compute | log(  (I3/I2)/t1  ) |
                array23 = zeros([length(x2)*length(x3), 1]);
                for pno = 1:length(x2)
                    for qno = 1:length(x3)
                        array23((pno-1)*length(x3)+qno,1) = abs( log( (inten(x3(qno))/inten(x2(pno)))/t1 ) );
                    end
                end
                [minlogR23, idx] = min(array23);

                % Decode selected pair indices (pno in x2, qno in x3)
                qno = mod(idx, length(x3));
                if 0==qno
                    qno = length(x3);
                end
                pno = (idx - qno)/length(x3) + 1;

                % 2.4 Accept the pair if:
                %   a) either M−1 is empty OR min(I in M−1) < I(M), AND
                %   b) the ratio deviation is within 'thr' (log-space check)
                if ( isempty(x1) || min(inten(x1)) < inten(x2(pno)) ) && minlogR23 <= log(thr)
                    % Store available channels:
                    if ~isempty(x1)
                        % M−1 kept as the *minimum* inside its window (original heuristic)
                        tmp_intens(1) = min(inten(x1));
                    end
                    tmp_intens(2) = inten(x2(pno));  % monoisotopic M
                    tmp_intens(3) = inten(x3(qno));  % M+1

                    % 2.5 Optionally add M+2 if ratio vs M+1 also within bounds (or forced)
                    if ~isempty(x4)
                        % Compare I(M+2)/I(M+1) to 't1' as a coarse check (original logic)
                        array34 = zeros([length(x4),1]);
                        for rno = 1:length(x4)
                            array34(rno,1) = abs( log( (inten(x4(rno))/inten(x3(qno)))/t1 ) );
                        end
                        [minlogR34, rno] = min(array34);

                        % Accept M+2 if (nC13==2) forced or ratio is within 'thr'
                        if 2==nC13 || minlogR34 <= log(thr)
                            tmp_intens(4) = inten(x4(rno));
                        end
                    end
                end
            end

            % 2.6 Save intensities for this scan (zeros if not accepted)
            c_ref_isointens(ino,1:nmz) = tmp_intens;
        end
    end
end
