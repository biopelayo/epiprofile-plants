function [cur_rts,cur_intens,cur_mono_isointens] = get_histone12(MS1_index,MS1_peaks,ptol,unitdiff,His,hno)
    % GET_HISTONE12 — Extracts RTs and intensities around a reference RT for a given histone peptide
    % INPUTS
    %   MS1_index  : [N x 2] index table for MS1 scans; column 2 must be retention times (RTs).
    %   MS1_peaks  : container/structure with MS1 peak lists compatible with GetProfiles (external).
    %   ptol       : precursor mass tolerance (same units used by GetProfiles; handled below).
    %   unitdiff   : nominal 13C isotopic mass difference (~1.003355... Da), used per charge for isotope m/z ladder.
    %   His        : structure holding histone peptide library and quant metadata (pep_mz, pep_ch, rt_ref, etc.).
    %   hno        : index of the peptide within His.* arrays to process.
    %
    % OUTPUTS
    %   cur_rts            : [1 x ncharge] RT positions selected for each charge state of the peptide.
    %   cur_intens         : [1 x ncharge] integrated areas (AUC) around the selected top apex for each charge.
    %   cur_mono_isointens : [T x 1] monoisotopic EIC intensities (from the first processed charge), used downstream.
    %
    % ASSUMPTIONS / CONTRACTS
    %   - His.pep_mz(hno,j) and His.pep_ch(hno,j) are defined for all charge states j used here.
    %   - His.rt_ref(hno) exists and is a scalar RT (reference neighborhood center).
    %   - GetProfiles, GetTopBottom, and get_area are available on the MATLAB path.
    %   - MS1_index(:,2) stores RTs sorted ascending; used to bracket a ±delta window.
    %
    % FAILURE MODES (not changed, only documented)
    %   - If the ±delta RT window has no points, "fallback" behavior below expands to whole chromatogram (Option 2).
    %   - A hard error path (Option 1) exists but is commented; you may enforce it for strict QC.

    [npep,ncharge] = size(His.pep_mz); %#ok
    % ^ Number of peptides (rows) and charge states (columns) defined in library.

    cur_rts = zeros([1,ncharge]);
    % ^ Will store the chosen RT per charge state.

    cur_intens = zeros([1,ncharge]);
    % ^ Will store integrated intensities per charge state.

    delta = 0.2;
    % ^ Half-window in minutes (±delta) around the reference RT. Narrow enough to focus on the target apex,
    %   but wide enough to account for minor RT drifts.

    %— ORIGINAL LINES (commented out by you above): they would fail if the window is empty
    % p  = find( MS1_index(:,2)>=His.rt_ref(hno)-delta );
    % rt_i1 = p(1);
    % pp = find( MS1_index(:,2)<=His.rt_ref(hno)+delta );
    % rt_i2 = pp(end);

    fprintf('DEBUG: hno=%d → RT_ref = %.3f\n', hno, His.rt_ref(hno));
    % ^ Debug print showing which peptide index we’re processing and its reference RT.

    %— NEW INDEXING: compute lower/upper indices bracketing the ±delta window safely
    lo_idx = find( MS1_index(:,2) >= His.rt_ref(hno) - delta );
    % ^ All scan indices whose RT is ≥ (RT_ref - delta); empty if none.

    hi_idx = find( MS1_index(:,2) <= His.rt_ref(hno) + delta );
    % ^ All scan indices whose RT is ≤ (RT_ref + delta); empty if none.

    %% OPCIÓN 1: Error controlado
    % Descomenta este bloque para que lance un error explícito si falta rango.
    % if isempty(lo_idx) || isempty(hi_idx)
    %     error('No hay puntos MS1 en ±%.2f de RT_ref=%.3f (hno=%d).', ...
    %           delta, His.rt_ref(hno), hno);
    % end
    %% que es cero o falta
    % Antes del cálculo de rt_i1 y rt_i2
    % fprintf('hno=%d: lo_idx = [%s], hi_idx = [%s]\n', ...
    %     hno, mat2str(lo_idx), mat2str(hi_idx));
    % 
    % rt_i1 = lo_idx(1);
    % rt_i2 = hi_idx(end);
    % if rt_i1 > rt_i2
    %     error('Ventana invertida: rt_i1 (%d) > rt_i2 (%d).', rt_i1, rt_i2);
    % end

    %% OPCIÓN 2: Fallback al cromatograma completo
    Descomenta este bloque para usar extremos si falta rango local.
    % ^ [WARNING] This line in your source is NOT commented and will cause a parse error if executed.
    %             I’m leaving it as-is because you asked for “tal cual”. If you want it commented, I’ll add '%'.

    if isempty(lo_idx)
        rt_i1 = 1;
        % ^ No lower bound found near RT_ref: fallback to beginning of MS1_index.
    else
        rt_i1 = lo_idx(1);
        % ^ First scan index meeting the lower RT bound.
    end
    if isempty(hi_idx)
        rt_i2 = size(MS1_index,1);
        % ^ No upper bound found: fallback to end of MS1_index.
    else
        rt_i2 = hi_idx(end);
        % ^ Last scan index meeting the upper RT bound.
    end
    % (Optional) sanity check: ensure start ≤ end
    if rt_i1 > rt_i2
        warning('rt_i1 (%d) > rt_i2 (%d): intercambiando.', rt_i1, rt_i2);
        tmp = rt_i1; rt_i1 = rt_i2; rt_i2 = tmp;
    end

    if ptol==100
        ptol = 10;
        % ^ Special case normalization: in some historical configs "100" meant "10 ppm" (or similar).
        %   Preserved behavior: if user passed 100, clamp to 10.
    end

    for jno=1:ncharge
        % === Loop over charge states defined for this peptide (columns of His.pep_mz) ===

        % -- Get MS1 monoisotopic profile for this charge state --
        c_mz = His.pep_mz(hno,jno);
        % ^ Theoretical precursor m/z for (hno, this charge).

        c_ch = His.pep_ch(hno,jno);
        % ^ Charge state value (z).

        % Build reference isotopic ladder around the monoisotope:
        %   [M-1], [M], [M+1], [M+2] at the correct spacing (unitdiff / z).
        c_ref_isomzs = [c_mz-unitdiff/c_ch c_mz c_mz+unitdiff/c_ch c_mz+2*unitdiff/c_ch];

        % Heuristic: allow one 13C shift for high-ptol and z≥3 (legacy behavior kept).
        if ptol>100 && c_ch>=3
            nC13 = 1;
        else
            nC13 = 0;
        end

        % Extract chromatographic profiles for the isotopic ladder within [rt_i1:rt_i2]
        % GetProfiles returns:
        %   c_isorts         : RT vector (sorted)
        %   c_ref_isointens  : [T x K] matrix of intensities for K isotopologues at each RT
        [c_isorts,c_ref_isointens] = GetProfiles(MS1_index,MS1_peaks,...
            c_ref_isomzs,c_ch,ptol,nC13,rt_i1:rt_i2);

        % We will use the monoisotopic trace (column j=2; indexing consistent with upstream code).
        j = 2;
        c_mono_isointens = c_ref_isointens(:,j);

        % Keep the first charge’s mono-isotopic trace for optional downstream reporting.
        if 1==jno
            cur_mono_isointens = c_mono_isointens;
        end

        % -- Peak picking: get apex and integration bounds from the mono trace --
        % GetTopBottom should return:
        %   nt        : indices of local maxima (apex candidates)
        %   nb        : integration boundaries (nb(k):left, nb(k)+1:right) aligned to nt(k)
        %   top1_idx  : best apex index among nt (by internal score)
        %   inten_sum : a per-peak area/intensity summary (vector aligned with nt)
        [nt,nb,top1_idx,inten_sum] = GetTopBottom(c_mono_isointens); %#ok

        % Restrict candidate apexes to those close to the library RT (±0.5 min window).
        flag = abs(c_isorts(nt) - His.rt_ref(hno)) <= 0.5;
        x = find(flag);

        if ~isempty(x)
            % If more than one candidate within ±0.5, choose the one with maximum area among them.
            [~,id] = max(inten_sum(x)); 
            top1_idx = x(id);

            % Convert from peak index to absolute position in c_isorts
            cur_pos = nt(top1_idx);

            % Store RT at chosen apex
            cur_rts(jno) = c_isorts(cur_pos);

            % If there are multiple candidates, rebuild the left/right integration bounds
            % using the outermost boundaries of the group to get a robust total area.
            if length(x)>=2
                nb = [nb(x(1)) nb(x(end)+1)];
            end

            % Integrate area around selected apex using get_area (consistent with EpiProfile)
            cur_intens(jno) = get_area(...
                c_isorts,c_ref_isointens,nb,cur_pos,c_mz,c_ch,...
                MS1_index,MS1_peaks,unitdiff,ptol);
        else
            % No acceptable apex within ±0.5 min of library RT: keep RT_ref and set area to 0.
            cur_rts(jno) = His.rt_ref(hno);
            cur_intens(jno) = 0;
        end
    end
end
