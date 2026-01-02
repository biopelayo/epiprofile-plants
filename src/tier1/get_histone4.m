function [cur_rts,cur_intens,cur_mono_isointens] = get_histone4( ...
    MS1_index, MS1_peaks, MS2_index, MS2_peaks, ptol, unitdiff, ...
    Mods, His, hno, special)
%% GET_HISTONE4
% Deconvolución en 4 componentes (p. ej., 4 isoformas PTM) a partir de:
%   - Área total MS1 (get_histone1)
%   - Ratios dirigidos por iones diagnósticos MS2
%
% Entradas (claves):
%   MS1_index(:,2): RT; MS1_index(:,3): offset a picos MS1
%   MS2_index(:,2): RT; MS2_index(:,4): precursor/centro; (:,6): instrumento; (:,7): offset a picos MS2
%   ptol: tolerancia ppm (si ==100 se fuerza a 10, por compatibilidad histórica)
%   unitdiff: masa 13C (~1.003355 Da)
%   His: struct con pep_mz, pep_ch, rt_ref, mod_short, outpath, outfile, ...
%   hno: índice de la PRIMERA isoforma del bloque (se usan hno:hno+3)
%   special.nDAmode: 1=DDA, 2=DIA (y otros flags como nsource/nsubtype)
%
% Salidas:
%   cur_rts            [4 x nCharge]   RT ápice por isoforma (replicado a todas las cargas)
%   cur_intens         [4 x nCharge]   Áreas por isoforma y carga (total MS1 * ratio_i)
%   cur_mono_isointens [num_MS1 x 4]   Sub-trazas mono separadas y alineadas a la rejilla MS1
%
% NOTA: Se mantiene la semántica original. Solo se mejora legibilidad y robustez.

% ---------------------------- Constantes/guardas ----------------------------
NPART = 4;                 % nº de isoformas a separar
DELTA = 0.5;               % margen temporal por defecto para ventana de trabajo
[num_MS1, ~] = size(MS1_index);
[npep, ncharge] = size(His.pep_mz); %#ok<ASGLU>

% Validación mínima del bloque hno:hno+3
if hno+3 > size(His.pep_mz,1) || hno < 1
    fprintf(1,'[get_histone4] Índice hno fuera de rango para 4 isoformas.\n');
    cur_rts = zeros(NPART, max(1,ncharge));
    cur_intens = zeros(NPART, max(1,ncharge));
    cur_mono_isointens = zeros(num_MS1, NPART);
    return;
end

% Inicialización de salidas
cur_rts            = zeros(NPART, ncharge);
cur_intens         = zeros(NPART, ncharge);
cur_mono_isointens = zeros(num_MS1, NPART);

% ---------------------- Cuantificación total (MS1) -------------------------
[h_rts, h_intens] = get_histone1(MS1_index, MS1_peaks, ptol, unitdiff, His, hno);
if h_intens(1) == 0
    % No hay señal MS1 utilizable → abortar (respeta el original)
    return;
end

% ---------------------- Estimación de ratios (MS2) -------------------------
[s_rts, ratio, cur_mono_isointens] = get_ratio_4iso( ...
    MS1_index, MS1_peaks, MS2_index, MS2_peaks, ptol, unitdiff, ...
    Mods, His, hno, h_rts, special);

% Replicar RTs por carga
cur_rts(1, :) = repmat(s_rts(1), 1, ncharge);
cur_rts(2, :) = repmat(s_rts(2), 1, ncharge);
cur_rts(3, :) = repmat(s_rts(3), 1, ncharge);
cur_rts(4, :) = repmat(s_rts(4), 1, ncharge);

% Repartir el área total por carga según ratios
cur_intens(1, :) = h_intens * ratio(1);
cur_intens(2, :) = h_intens * ratio(2);
cur_intens(3, :) = h_intens * ratio(3);
cur_intens(4, :) = h_intens * ratio(4);

% ----------------------- Figura de control (PDF) ---------------------------
try
    set(gcf, 'visible', 'off');

    % Traza total y sub-trazas (sobre la ventana utilizada)
    c_mz  = His.pep_mz(hno,1);
    c_ch  = His.pep_ch(hno,1);
    if ptol == 100, ptol = 10; end
    nC13  = (ptol > 100 && c_ch >= 3);
    rt_i1 = find(MS1_index(:,2) >= His.rt_ref(hno)   - DELTA, 1, 'first');
    rt_i2 = find(MS1_index(:,2) <= His.rt_ref(hno+3) + DELTA, 1, 'last');
    c_ref_isomzs = [c_mz - unitdiff/c_ch, c_mz, c_mz + unitdiff/c_ch, c_mz + 2*unitdiff/c_ch];
    [c_isorts, c_ref_isointens] = GetProfiles(MS1_index, MS1_peaks, c_ref_isomzs, c_ch, ptol, nC13, rt_i1:rt_i2);
    mono = c_ref_isointens(:,2);                   % traza mono experimental en ventana
    n_rt = c_isorts(rt_i1:rt_i2);

    % Extraer las cuatro sub-trazas desde la matriz de salida
    n_inten1 = mono(rt_i1:rt_i2);                  % total experimental (referencia visual)
    n_inten2 = cur_mono_isointens(rt_i1:rt_i2,1);  % iso 1
    n_inten3 = cur_mono_isointens(rt_i1:rt_i2,2);  % iso 2
    n_inten4 = cur_mono_isointens(rt_i1:rt_i2,3);  % iso 3
    n_inten5 = cur_mono_isointens(rt_i1:rt_i2,4);  % iso 4

    plot(n_rt, n_inten1, 'linestyle','-',  'linewidth',2, 'color','k'); hold on;
    plot(n_rt, n_inten2, 'linestyle','--', 'linewidth',2, 'color','r');
    plot(n_rt, n_inten3, 'linestyle','-.', 'linewidth',2, 'color','b');
    plot(n_rt, n_inten4, 'linestyle','--', 'linewidth',2, 'color','g');
    plot(n_rt, n_inten5, 'linestyle','-.', 'linewidth',2, 'color','m');
    xlabel('time (min)'); ylabel('intensity');

    legend('experiment', His.mod_short{hno}, His.mod_short{hno+1}, His.mod_short{hno+2}, His.mod_short{hno+3});

    r1 = floor(100 * ratio(1));
    r2 = floor(100 * ratio(2));
    r3 = floor(100 * ratio(3));
    r4 = max(0, 100 - r1 - r2 - r3);
    title(sprintf('%s/%s/%s/%s:%d%%:%d%%:%d%%:%d%%', ...
        His.mod_short{hno}, His.mod_short{hno+1}, His.mod_short{hno+2}, His.mod_short{hno+3}, ...
        r1, r2, r3, r4));

    xx = strfind(His.mod_short{hno}, 'ac');           % nº de "ac" en el primer label (como en el original)
    out_file1 = fullfile(His.outpath, ...
        ['Iso_', His.outfile, '_', num2str(length(xx)), 'ac', '_', His.mod_short{hno}, '.pdf']);
    print('-dpdf', out_file1);
    close();
catch ME
    % Nunca abortar por un fallo de figura
    fprintf(1,'[get_histone4] Aviso al generar PDF: %s\n', ME.message);
    close all force;
end

% =========================================================================
% FUNCIONES ANIDADAS
% =========================================================================
function [s_rts, ratio, cur_mono_isointens] = get_ratio_4iso( ...
    MS1_index, MS1_peaks, MS2_index, MS2_peaks, ptol, unitdiff, ...
    Mods, His, hno, h_rts, special) %#ok<INUSD>
% Calcula los RT de ápice de cada isoforma y sus ratios globales,
% además de construir las sub-trazas mono por isoforma.

    npart  = NPART;
    s_rts  = zeros(1, npart);
    ratio  = zeros(1, npart);
    num_MS1 = size(MS1_index,1);
    cur_mono_isointens = zeros(num_MS1, npart);

    % Ventana temporal basada en RTs de referencia de la 1ª y 4ª isoformas
    pL = find(MS1_index(:,2) >= His.rt_ref(hno)   - DELTA, 1, 'first');
    pR = find(MS1_index(:,2) <= His.rt_ref(hno+3) + DELTA, 1, 'last');

    % Traza mono experimental en la ventana
    c_mz = His.pep_mz(hno,1);
    c_ch = His.pep_ch(hno,1);
    if ptol == 100, ptol = 10; end
    nC13 = (ptol > 100 && c_ch >= 3);
    seeds = [c_mz - unitdiff/c_ch, c_mz, c_mz + unitdiff/c_ch, c_mz + 2*unitdiff/c_ch];
    [c_isorts, c_ref_isointens] = GetProfiles(MS1_index, MS1_peaks, seeds, c_ch, ptol, nC13, pL:pR);
    mono = c_ref_isointens(:,2);
    IX   = pL:pR;
    rt1  = c_isorts(IX(1));      % ventana efectiva en RT (puede refinarse en DIA)
    rt2  = c_isorts(IX(end));

    % En DIA, restringir a la ventana cuyo centro coincide con c_mz
    if special.nDAmode == 2
        premzs  = unique(MS2_index(:,4));
        [~, ii] = min(abs(premzs - c_mz));
        target  = premzs(ii);

        flag = false(size(MS2_index,1),1);
        i1   = find(MS2_index(:,2) >= rt1, 1, 'first');
        i2   = find(MS2_index(:,2) <= rt2,  1, 'last');
        if isempty(i1) || isempty(i2), return; end

        for i = i1:i2
            flag(i) = (MS2_index(i,4) == target);
        end

        ms2pos = find(flag);
        if isempty(ms2pos), return; end

        % Remapear a índices MS1 para fijar rt1/rt2 al rango real observado en DIA
        IXdia = zeros(numel(ms2pos),1);
        for k = 1:numel(ms2pos)
            IXdia(k) = find(MS1_index(:,1) == MS2_index(ms2pos(k),1), 1, 'first');
        end
        rt1 = c_isorts(IXdia(1));
        rt2 = c_isorts(IXdia(end));
    end

    % Emparejar MS2 y obtener proporciones por escaneo
    [ms2pos, ms2ratios] = MatchMS2(MS2_index, MS2_peaks, c_mz, c_ch, ptol, unitdiff, ...
                                   Mods, His, hno, rt1, rt2, special);

    % Si hay muy poca evidencia MS2, fallback uniforme
    nlen = numel(IX);
    if numel(ms2pos) <= 3
        if nlen <= 3, return; end
        s_rts(1:npart) = [His.rt_ref(hno), His.rt_ref(hno+1), His.rt_ref(hno+2), His.rt_ref(hno+3)];
        ratio(1:npart) = repmat(1/npart, 1, npart);

        n_inten = mono(IX) / npart;
        if special.nDAmode ~= 2
            cur_mono_isointens(:,1) = [zeros(pL-1,1); n_inten; zeros(num_MS1-pR,1)];
            cur_mono_isointens(:,2) = cur_mono_isointens(:,1);
            cur_mono_isointens(:,3) = cur_mono_isointens(:,1);
            cur_mono_isointens(:,4) = cur_mono_isointens(:,1);
        else
            xx = c_isorts(IX(1):IX(end));
            yy = spline(c_isorts(IX), n_inten, xx);
            pad = [zeros(pL-1,1); yy; zeros(num_MS1-pR,1)];
            cur_mono_isointens = [pad, pad, pad, pad];
        end
        return;
    end

    % Construir matriz de ratios por escaneo (nlen x 4)
    rmat   = zeros(nlen, npart);
    scans2 = MS2_index(ms2pos,1);

    for i = 1:nlen
        scan = MS1_index(IX(i),1);
        [tf, loc] = ismember(scan, scans2);
        if tf
            rmat(i,1:npart) = ms2ratios(1:npart, loc).';
        end
    end

    % Suavizado y normalización fila a fila
    for k = 1:npart
        rmat(:,k) = smooth(rmat(:,k), 3);
    end
    rowSum = sum(rmat,2);
    rowSum(rowSum == 0) = eps;
    for k = 1:npart
        rmat(:,k) = rmat(:,k) ./ rowSum;
    end

    % Separación de la traza mono experimental en 4 sub-trazas
    n_rt    = c_isorts(IX);
    mono_v  = mono(IX);
    i2 = mono_v .* rmat(:,1);
    i3 = mono_v .* rmat(:,2);
    i4 = mono_v .* rmat(:,3);
    i5 = mono_v .* rmat(:,4);

    % Ápices por isoforma
    [~, x1] = max(i2);
    [~, x2] = max(i3);
    [~, x3] = max(i4);
    [~, x4] = max(i5);

    % Áreas por *spline* (paso fijo)
    xx = n_rt(1):0.005:n_rt(end);
    a1 = sum(abs(spline(n_rt, i2, xx)));
    a2 = sum(abs(spline(n_rt, i3, xx)));
    a3 = sum(abs(spline(n_rt, i4, xx)));
    a4 = sum(abs(spline(n_rt, i5, xx)));
    S  = eps + a1 + a2 + a3 + a4;

    s_rts(1:npart) = [n_rt(x1), n_rt(x2), n_rt(x3), n_rt(x4)];
    ratio(1:npart) = [a1, a2, a3, a4] ./ S;

    % Volcado de sub-trazas a la rejilla global MS1
    if special.nDAmode ~= 2
        cur_mono_isointens(:,1) = [zeros(pL-1,1); i2; zeros(num_MS1-pR,1)];
        cur_mono_isointens(:,2) = [zeros(pL-1,1); i3; zeros(num_MS1-pR,1)];
        cur_mono_isointens(:,3) = [zeros(pL-1,1); i4; zeros(num_MS1-pR,1)];
        cur_mono_isointens(:,4) = [zeros(pL-1,1); i5; zeros(num_MS1-pR,1)];
    else
        % En DIA, re-muestrear sobre la malla original (suaviza bordes)
        xxg = c_isorts(IX(1):IX(end));
        y1  = spline(n_rt, i2, xxg);
        y2  = spline(n_rt, i3, xxg);
        y3  = spline(n_rt, i4, xxg);
        y4  = spline(n_rt, i5, xxg);
        cur_mono_isointens(:,1) = [zeros(pL-1,1); y1; zeros(num_MS1-pR,1)];
        cur_mono_isointens(:,2) = [zeros(pL-1,1); y2; zeros(num_MS1-pR,1)];
        cur_mono_isointens(:,3) = [zeros(pL-1,1); y3; zeros(num_MS1-pR,1)];
        cur_mono_isointens(:,4) = [zeros(pL-1,1); y4; zeros(num_MS1-pR,1)];
    end
end

function [ms2pos, ms2ratios] = MatchMS2( ...
    MS2_index, MS2_peaks, c_mz, c_ch, ptol, unitdiff, ...
    Mods, His, hno, rt1, rt2, special)
% Empareja escaneos MS2 relevantes y calcula proporciones por isoforma
% usando iones diagnósticos (tres comparaciones secuenciales 1vs2, 2vs3, 3vs4).

    % Señal NH-mass según origen/subtipo
    nhmass = (special.nsource == 4) && (special.nsubtype ~= 0) && (special.nsubtype ~= 2);

    % --- Selección de MS2 relevantes (DDA/DIA) ---
    num_MS2 = size(MS2_index,1);
    if special.nDAmode == 1
        % DDA: aceptar precursores M, M+1, M+2 (en m/z) dentro de ptol ppm
        sets = [0 1 2];
        mzs  = c_mz + sets * unitdiff / c_ch;
        flag = false(num_MS2,1);

        i1 = find(MS2_index(:,2) >= rt1, 1, 'first');
        i2 = find(MS2_index(:,2) <= rt2,  1, 'last');
        if isempty(i1) || isempty(i2)
            ms2pos = []; ms2ratios = []; return;
        end

        for i = i1:i2
            cen_mz = MS2_index(i,4);
            ok = any(abs(mzs - cen_mz) < ptol * cen_mz * 1e-6);
            flag(i) = ok;
        end
        ms2pos = find(flag);
        if isempty(ms2pos), ms2ratios = []; return; end

    elseif special.nDAmode == 2
        % DIA: ventana cuyo centro coincide con c_mz
        premzs = unique(MS2_index(:,4));
        [~, ii] = min(abs(premzs - c_mz));
        target  = premzs(ii);

        flag = false(num_MS2,1);
        i1   = find(MS2_index(:,2) >= rt1, 1, 'first');
        i2   = find(MS2_index(:,2) <= rt2,  1, 'last');
        if isempty(i1) || isempty(i2)
            ms2pos = []; ms2ratios = []; return;
        end

        for i = i1:i2
            flag(i) = (MS2_index(i,4) == target);
        end
        ms2pos = find(flag);
        if isempty(ms2pos), ms2ratios = []; return; end

    else
        ms2pos = []; ms2ratios = []; return;
    end

    % --- Tolerancia y tipo de fragmentación por instrumento ---
    instruments = MS2_index(ms2pos,6); % {CIDIT,CIDFT,ETDIT,ETDFT,HCDIT,HCDFT}
    singleInstr = (numel(unique(instruments)) == 1);
    if singleInstr
        inst = instruments(1);
        ActiveType = ternary( any(inst == [3 4]), 'ETD', 'CID' );
        tol = ternary( mod(inst,2)==1, 0.4, 0.02 ); % IT vs FT
        % Listas de iones diagnósticos:
        if nhmass
            [K11,K12] = get_key_ionsH(His,hno,  hno+1, Mods, ActiveType);
            [K21,K22] = get_key_ionsH(His,hno+1,hno+2, Mods, ActiveType);
            [K31,K32] = get_key_ionsH(His,hno+2,hno+3, Mods, ActiveType);
        else
            [K11,K12] = get_key_ions (His,hno,  hno+1, Mods, ActiveType);
            [K21,K22] = get_key_ions (His,hno+1,hno+2, Mods, ActiveType);
            [K31,K32] = get_key_ions (His,hno+2,hno+3, Mods, ActiveType);
        end
    end

    % --- Ratios por escaneo MS2 seleccionado ---
    index = [1; MS2_index(1:num_MS2,7)];
    ms2ratios = zeros(NPART, numel(ms2pos)); % filas: isoformas 1..4

    for i = 1:numel(ms2pos)
        cno = ms2pos(i);
        for pno = cno % (en el original se llegó a probar ±1; aquí se mantiene el central)
            if pno < 1 || pno > num_MS2, continue; end

            % Si hay mezcla de instrumentos, recalcular tol y ActiveType por escaneo
            if ~singleInstr
                inst = MS2_index(pno,6);
                ActiveType = ternary( any(inst == [3 4]), 'ETD', 'CID' );
                tol = ternary( mod(inst,2)==1, 0.4, 0.02 );
                if nhmass
                    [K11,K12] = get_key_ionsH(His,hno,  hno+1, Mods, ActiveType);
                    [K21,K22] = get_key_ionsH(His,hno+1,hno+2, Mods, ActiveType);
                    [K31,K32] = get_key_ionsH(His,hno+2,hno+3, Mods, ActiveType);
                else
                    [K11,K12] = get_key_ions (His,hno,  hno+1, Mods, ActiveType);
                    [K21,K22] = get_key_ions (His,hno+1,hno+2, Mods, ActiveType);
                    [K31,K32] = get_key_ions (His,hno+2,hno+3, Mods, ActiveType);
                end
            end

            IX2 = index(pno):index(pno+1)-1;
            mz  = MS2_peaks(IX2,1);
            in  = MS2_peaks(IX2,2);

            % Tres comparaciones secuenciales → a, ab, abc → reconstrucción de [a b c d]
            a   = get_ratio(mz, in, tol, K11, K12);
            ab  = get_ratio(mz, in, tol, K21, K22);
            abc = get_ratio(mz, in, tol, K31, K32);

            d = 1 - abc;
            b = max(ab - a, 0);
            c = 1 - a - b - d;

            if a + d > 1
                d = 1 - a; b = 0; c = 0;
            elseif a + d + b > 1
                b = 1 - a - d; c = 0;
            end

            X = [a; b; c; d];
            % Mantener preferencia por vectores menos degenerados (menos ceros)
            if nnz(ms2ratios(:,i) == 0) > nnz(X == 0)
                ms2ratios(:,i) = X;
            end
        end
    end

    % Si todos los escaneos dieron d==1 (degen.), normalizar a [0 0 0.01 0.99] en esos puntos
    x = find(ms2ratios(NPART,:) < 1); %#ok<FNDSB>
    if isempty(x)
        ms2pos    = [];
        ms2ratios = [];
    else
        xx = find(ms2ratios(NPART,:) == 1);
        for ino = 1:numel(xx)
            ms2ratios(:,xx(ino)) = [0; 0; 0.01; 0.99];
        end
    end
end

function ratio1 = get_ratio(mz, inten, tol, K1, K2)
% Suma de intensidades en ventanas ±tol para dos listas diagnósticas K1 y K2
    i1 = zeros(1, numel(K1));
    i2 = zeros(1, numel(K2));
    for j = 1:numel(K1)
        ix1 = find(abs(mz - K1(j)) <= tol);
        ix2 = find(abs(mz - K2(j)) <= tol);
        if ~isempty(ix1)
            [~, x1] = min(abs(mz(ix1) - K1(j)));
            i1(j) = inten(ix1(x1));
        end
        if ~isempty(ix2)
            [~, x2] = min(abs(mz(ix2) - K2(j)));
            i2(j) = inten(ix2(x2));
        end
    end
    % Proporción agregada
    ratio1 = sum(i1) / (eps + sum(i1) + sum(i2));
end

function out = ternary(cond, A, B)
% Operador ternario simple
    if cond, out = A; else, out = B; end
end
end

