function [cur_rts,cur_intens,cur_mono_isointens] = get_histone12(MS1_index,MS1_peaks,ptol,unitdiff,His,hno)
% get_histone12 — Robust extraction of monoisotopic MS1 traces in a tight RT window (±0.2 min),
%                 with defensive indexing and optional fallback to the full chromatogram range.
% -------------------------------------------------------------------------------------------------
% INPUTS
%   MS1_index   : [nScans x 2] matrix; col 2 contains retention times (minutes).
%   MS1_peaks   : project-specific structure/array used by GetProfiles to read MS1 intensities.
%   ptol        : mass tolerance parameter; if exactly 100, it is internally remapped to 10.
%   unitdiff    : isotopic spacing per charge (≈ 1.00335 u divided by z).
%   His         : struct with histone peptide metadata:
%                   - pep_mz(hno, j) : theoretical m/z for peptide hno at charge state j
%                   - pep_ch(hno, j) : charge state for peptide hno at column j
%                   - rt_ref(hno)    : reference RT (minutes) for peptide hno
%   hno         : integer index of the peptide row within His.pep_mz/pep_ch and His.rt_ref.
%
% OUTPUTS
%   cur_rts            : [1 x ncharge] selected apex RT per charge state.
%   cur_intens         : [1 x ncharge] integrated area per charge state.
%   cur_mono_isointens : monoisotopic intensity time series for the FIRST charge only (diagnostics).
%
% NOTES
%   - Compared to the minimal variant, this version prints DEBUG information and guards against
%     empty index ranges by *falling back* to the full chromatogram (start..end) when needed.
%   - Peak finding uses GetTopBottom (not GetTopBottom11). Candidate filtering still uses ±0.5 min
%     around His.rt_ref even though the extraction window is ±0.2.
% -------------------------------------------------------------------------------------------------

    % Determine number of peptides (rows) and number of charges (columns).
    [npep,ncharge] = size(His.pep_mz); %#ok

    % Preallocate outputs (one value per charge).
    cur_rts    = zeros([1,ncharge]);
    cur_intens = zeros([1,ncharge]);

    % Narrow half-window (minutes) for local RT extraction.
    delta = 0.2;

    % ---- ORIGINAL LINES (commented in this robust version) -----------------
    % p  = find( MS1_index(:,2)>=His.rt_ref(hno)-delta );
    % rt_i1 = p(1);
    % pp = find( MS1_index(:,2)<=His.rt_ref(hno)+delta );
    % rt_i2 = pp(end);
    % -----------------------------------------------------------------------

    % Emit DEBUG header for this peptide index and its reference RT.
    fprintf('DEBUG: hno=%d → RT_ref = %.3f\n', hno, His.rt_ref(hno));

    % Compute candidate scan indices on the lower/upper side of the RT window.
    lo_idx = find( MS1_index(:,2) >= His.rt_ref(hno) - delta );
    fprintf('DEBUG: lo_idx tiene %d elementos:\n', numel(lo_idx));
    disp(lo_idx.');  % print as row for readability

    hi_idx = find( MS1_index(:,2) <= His.rt_ref(hno) + delta );
    fprintf('DEBUG: hi_idx tiene %d elementos:\n', numel(hi_idx));
    disp(hi_idx.');

    % OPTIONAL strict error path (kept commented): enforce local RT coverage.
    % if isempty(lo_idx) || isempty(hi_idx)
    %     error('No hay puntos MS1 en ±%.2f de RT_ref=%.3f (hno=%d).', ...
    %           delta, His.rt_ref(hno), hno);
    % end

    % OPTIONAL debug of raw indices (kept commented).
    % fprintf('hno=%d: lo_idx = [%s], hi_idx = [%s]\n', ...
    %     hno, mat2str(lo_idx), mat2str(hi_idx));
    %
    % rt_i1 = lo_idx(1);
    % rt_i2 = hi_idx(end);
    % if rt_i1 > rt_i2
    %     error('Ventana invertida: rt_i1 (%d) > rt_i2 (%d).', rt_i1, rt_i2);
    % end

    % ---- FALLBACK STRATEGY: use full chromatogram bounds if local window is empty.
    if ~isempty(lo_idx) && ~isempty(hi_idx)
        fprintf('DEBUG: rt_i1 = %d (lo_idx(1)), rt_i2 = %d (hi_idx(end))\n', ...
            lo_idx(1), hi_idx(end));
    end

    % If lower-bound set is empty, default to the first scan; otherwise pick first candidate.
    if isempty(lo_idx)
        rt_i1 = 1;
    else
        rt_i1 = lo_idx(1);
    end

    % If upper-bound set is empty, default to the last scan; otherwise pick last candidate.
    if isempty(hi_idx)
        rt_i2 = size(MS1_index,1);
    else
        rt_i2 = hi_idx(end);
    end

    % Ensure the indices define a forward window; if inverted, swap and warn.
    if rt_i1 > rt_i2
        warning('rt_i1 (%d) > rt_i2 (%d): intercambiando.', rt_i1, rt_i2);
        tmp = rt_i1; rt_i1 = rt_i2; rt_i2 = tmp;
    end

    % Legacy tolerance remapping: 100 means "use 10" internally (stricter extraction).
    if ptol==100
        ptol = 10;
    end

    % Iterate through all charge states available for this peptide.
    for jno=1,ncharge

        % --- Build the four isotopologue anchors: [M-1, M, M+1, M+2] ---
        c_mz = His.pep_mz(hno,jno);  % theoretical monoisotopic m/z
        c_ch = His.pep_ch(hno,jno);  % charge state (z)
        c_ref_isomzs = [ ...
            c_mz - unitdiff/c_ch, ...  % M-1 (symmetric reference)
            c_mz, ...                  % M (monoisotopic)
            c_mz + unitdiff/c_ch, ...  % M+1
            c_mz + 2*unitdiff/c_ch ... % M+2
        ];

        % Heuristic: coarse tolerance & high charge may require different envelope handling.
        if ptol>100 && c_ch>=3
            nC13 = 1;
        else
            nC13 = 0;
        end

        % --- Extract aligned profiles within the chosen scan range [rt_i1:rt_i2] ---
        [c_isorts,c_ref_isointens] = GetProfiles( ...
            MS1_index,MS1_peaks, ...
            c_ref_isomzs,c_ch,ptol,nC13, ...
            rt_i1:rt_i2);

        % The monoisotopic intensity trace is the second column by convention.
        j = 2;
        c_mono_isointens = c_ref_isointens(:,j);

        % Keep the full monoisotopic trace only for the first charge (diagnostics).
        if 1==jno
            cur_mono_isointens = c_mono_isointens;
        end

        % --- Peak candidates and integration boundaries via GetTopBottom (variant) ---
        [nt,nb,top1_idx,inten_sum] = GetTopBottom(c_mono_isointens); %#ok

        % Candidate filter: keep tops within ±0.5 min of the reference RT.
        flag = abs(c_isorts(nt) - His.rt_ref(hno)) <= 0.5;
        x = find(flag);

        if ~isempty(x)
            % Among nearby candidates, choose the one with largest integrated proxy (area).
            [~,id]   = max(inten_sum(x));    %#ok
            top1_idx = x(id);                % index into nt of the winning candidate
            cur_pos  = nt(top1_idx);         % position into c_isorts (time grid)

            % Report apex RT for this charge state.
            cur_rts(jno) = c_isorts(cur_pos);

            % If multiple candidates, coalesce integration window from first to last.
            if length(x)>=2
                nb = [nb(x(1)) nb(x(end)+1)];  % NOTE: assumes nb has a trailing boundary.
            end

            % Integrate area for this charge using the selected bounds and apex.
            cur_intens(jno) = get_area( ...
                c_isorts,c_ref_isointens,nb,cur_pos, ...
                c_mz,c_ch,MS1_index,MS1_peaks,unitdiff,ptol);

        else
            % No valid candidate near rt_ref: keep RT at reference and area = 0.
            cur_rts(jno)    = His.rt_ref(hno);
            cur_intens(jno) = 0;
        end
    end
end
