function [rts, top1_rt, inten_sum, top1_inten_sum] = get_rts(MS1_index, MS1_peaks, ptol, unitdiff, His, hno, nsplit, t1, t2)
%%
% GET_RTS Extrae tiempos de retención (RT) e intensidades para un péptido/histona
% en el rango [t1, t2], usando siempre la primera carga (His.pep_ch(hno,1)).  
%
%   [rts, top1_rt, inten_sum, top1_inten_sum] = get_rts( ...
%       MS1_index, MS1_peaks, ptol, unitdiff, His, hno, nsplit, t1, t2 )
%
%   ENTRADAS:
%     • MS1_index      : Matriz N×M donde la segunda columna (MS1_index(:,2)) contiene
%                        los tiempos de retención de cada escaneo MS1 (en minutos), en 
%                        orden ascendente.
%     • MS1_peaks      : Matriz de intensidades MS1 asociadas a cada fila de MS1_index.
%     • ptol           : Tolerancia en m/z para la extracción de perfiles. Si vale 100,
%                        se ajusta internamente a 10 (caso “por defecto”).
%     • unitdiff       : Separación isotópica en masa (aprox. 1 Da dividido entre carga).
%     • His            : Estructura con datos del péptido/histona:
%         – His.pep_mz  : Matriz [npep × ncharge] de valores teóricos de m/z por carga.
%         – His.pep_ch  : Matriz [npep × ncharge] con número de carga de cada péptido.
%         – His.rt_ref  : Vector [npep × 1] con RT de referencia (no se usa directamente aquí,
%                        pero se documenta para entender el flujo general de procesado).
%     • hno            : Índice (posición) del péptido/histona en His a procesar.
%     • nsplit         : Si es 1, se llama a GetTopBottom; si >1, se llama a GetTopBottom11.
%     • t1, t2         : Límite inferior y superior (en minutos) del rango de RT a procesar.
%
%   SALIDAS:
%     • rts             : Vector con todos los RT detectados (c_isorts(nt)).
%     • top1_rt         : RT correspondiente al pico principal (c_isorts(nt(top1_idx))).
%     • inten_sum       : Vector con la suma de intensidades alrededor de cada pico candidato.
%     • top1_inten_sum  : Suma de intensidades correspondiente a top1_rt.
%
%   NOTA SOBRE “Array indices must be positive integers or logical values”:
%     Este error aparece cuando se intenta indexar un arreglo con un índice que:
%       - Es 0 o negativo.
%       - Es un decimal (no entero).
%       - Excede el tamaño del arreglo.
%       - El vector de índices está vacío y se intenta usar p(1) o pp(end).
%
%     En particular, en este bloque:
%       p  = find(MS1_index(:,2) >= t1);    % busca RT ≥ t1
%       rt_i1 = p(1);                       % falla si p == []
%       pp = find(MS1_index(:,2) <= t2);    % busca RT ≤ t2
%       rt_i2 = pp(end);                    % falla si pp == []
%
%     Para evitar dicho error y, en este caso, “ampliar un poco la ventana cuando t2
%     quede por debajo del mínimo RT disponible”, añadimos estos chequeos:
%       • Si t2 < min_RT_MS1, ajustamos t2 = min_RT_MS1 (con warning).
%       • Si t1 < min_RT_MS1, ajustamos t1 = min_RT_MS1 (con warning).
%
%   FLUJO DE LA FUNCIÓN:
%     1) Validar que [t1, t2] solape con los RT disponibles en MS1_index(:,2), con posible
%        “ampliación” si t2 < RT mínimo.
%     2) Extraer m/z e isótopos de His.pep_mz(hno,1) y su carga His.pep_ch(hno,1).
%     3) Localizar índices rt_i1 y rt_i2 para t1 ≤ RT ≤ t2, validando p y pp.
%     4) Llamar a GetProfiles para obtener c_isorts y c_ref_isointens en ese rango.
%     5) Seleccionar la columna de intensidades monoisotópicas.
%     6) Llamar a GetTopBottom o GetTopBottom11 según nsplit, obteniendo nt, nb, top1_idx, inten_sum.
%     7) Validar que nt no esté vacío y que nt contenga índices válidos para c_isorts.
%     8) Construir rts = c_isorts(nt), top1_rt = c_isorts(nt(top1_idx)), top1_inten_sum = inten_sum(top1_idx).
%     9) Si top1_rt < 4, devolver todo vacío (criterio de filtrado).
%
%   EJEMPLO DE USO:
%     % Supongamos que 'His' ya contiene pep_mz, pep_ch, rt_ref, etc.
%     valid_hno = find(His.rt_ref > 0);  % Solo procesar péptidos con rt_ref > 0
%     for k = 1:numel(valid_hno)
%         hno = valid_hno(k);
%         t1 = His.rt_ref(hno) - 0.5;
%         t2 = His.rt_ref(hno) + 0.5;
%         [rts, top1_rt, inten_sum, top1_inten_sum] = get_rts( ...
%             MS1_index, MS1_peaks, ptol, unitdiff, His, hno, 2, t1, t2 );
%         % … almacenar o procesar resultados …
%     end
%

    %% 1) VALIDAR Y, SI ES NECESARIO, AMPLIAR [t1, t2] PARA CUBRIR RT MÍNIMO DE MS1
    % Si t2 < RT mínimo de MS1, ajustamos t2 hacia abajo para no quedar “fuera de rango”.
    % Similarmente, si t1 < RT mínimo, lo ajustamos al mismo valor mínimo.
    % De este modo evitamos el error “ningún RT_MS1 <= t2” cuando t2 sea demasiado pequeño.
    min_RT_MS1 = MS1_index(1, 2);        % El valor más bajo de RT en MS1_index
    max_RT_MS1 = MS1_index(end, 2);      % El valor más alto de RT en MS1_index

    % 1.1) Si t1 > RT máximo, no hay solapamiento: retornamos vacío
    if t1 > max_RT_MS1
        rts            = [];
        top1_rt        = [];
        inten_sum      = [];
        top1_inten_sum = [];
        return;
    end

    % 1.2) Si t2 < 0 o t2 < t1, no es un rango válido → retornar vacío
    if t2 < 0 || t2 < t1
        rts            = [];
        top1_rt        = [];
        inten_sum      = [];
        top1_inten_sum = [];
        return;
    end

    % 1.3) Si t2 está por debajo del RT mínimo, ampliamos t2 = min_RT_MS1 (con warning)
    if t2 < min_RT_MS1
        warning( ...
            'get_rts: t2 (%.4g) < RT mínimo MS1 (%.4g). Ajustando t2 = %.4g.', ...
            t2, min_RT_MS1, min_RT_MS1 ...
        );
        t2 = min_RT_MS1;
    end

    % 1.4) Si t1 está por debajo del RT mínimo, lo igualamos a min_RT_MS1 (con warning)
    if t1 < min_RT_MS1
        warning( ...
            'get_rts: t1 (%.4g) < RT mínimo MS1 (%.4g). Ajustando t1 = %.4g.', ...
            t1, min_RT_MS1, min_RT_MS1 ...
        );
        t1 = min_RT_MS1;
    end

    %% 2) EXTRAER m/z TEÓRICO E ISÓTOPOS PARA LA CARGA 1
    c_mz = His.pep_mz(hno, 1);
    c_ch = His.pep_ch(hno, 1);

    % Vector de m/z de isótopos: [-1, 0, +1, +2] unidades de masa normalizadas por carga
    c_ref_isomzs = [ ...
        c_mz - unitdiff / c_ch, ...
        c_mz, ...
        c_mz + unitdiff / c_ch, ...
        c_mz + 2 * unitdiff / c_ch ...
    ];

    % Si ptol == 100, se reasigna internamente a 10
    if ptol == 100
        ptol = 10;
    end

    % Determinar si se considera isótopo +1 (nC13 = 1) o no (nC13 = 0)
    if ptol > 100 && c_ch >= 3
        nC13 = 1;
    else
        nC13 = 0;
    end

    %% 3) LOCALIZAR ÍNDICES rt_i1 Y rt_i2 EN MS1_index(:,2) PARA [t1, t2]
    % 3.1) Buscar p = índices donde RT ≥ t1
    p = find(MS1_index(:, 2) >= t1);
    if isempty(p)
        % Si no hay ningún escaneo con RT ≥ t1 (incluso tras ajuste), retornamos vacío
        rts            = [];
        top1_rt        = [];
        inten_sum      = [];
        top1_inten_sum = [];
        return;
    end
    rt_i1 = p(1);  % Primer índice válido

    % 3.2) Buscar pp = índices donde RT ≤ t2
    pp = find(MS1_index(:, 2) <= t2);
    if isempty(pp)
        % Teóricamente no debería ocurrir luego del ajuste, pero lo chequeamos:
        % Si aún así no hay escaneo RT ≤ t2, retornamos vacío
        rts            = [];
        top1_rt        = [];
        inten_sum      = [];
        top1_inten_sum = [];
        return;
    end
    rt_i2 = pp(end);  % Último índice válido

    %% 4) EXTRAER PERFIL MS1 ENTRE rt_i1 Y rt_i2
    [c_isorts, c_ref_isointens] = GetProfiles( ...
        MS1_index, MS1_peaks, c_ref_isomzs, c_ch, ptol, nC13, rt_i1:rt_i2 ...
    );

    % Tomar solo la columna 2 (monoisótopo)
    c_mono_isointens = c_ref_isointens(:, 2);

    %% 5) IDENTIFICAR “TOPS” Y “BOTTOMS” EN EL PERFIL MONOISOTÓPICO
    if nsplit == 1
        [nt, nb, top1_idx, inten_sum] = GetTopBottom(c_mono_isointens); %#ok
    else
        [nt, nb, top1_idx, inten_sum] = GetTopBottom11(c_mono_isointens); %#ok
    end

    % 5.1) Si nt está vacío, no se detectaron picos “top”: retornamos vacío
    if isempty(nt)
        rts            = [];
        top1_rt        = [];
        inten_sum      = [];
        top1_inten_sum = [];
        return;
    end

    % 5.2) Comprobar que todos los valores de nt sean índices válidos para c_isorts
    if any(nt < 1) || any(nt > numel(c_isorts)) || any(floor(nt) ~= nt)
        error('get_rts: nt contiene índices inválidos (fuera de [1, %d] o no enteros).', numel(c_isorts));
    end

    % 5.3) Construir rts: RT reales de cada pico detectado (c_isorts(nt))
    rts = c_isorts(nt);

    % 5.4) Verificar que top1_idx sea un entero dentro de [1, numel(nt)]
    if top1_idx < 1 || top1_idx > numel(nt) || floor(top1_idx) ~= top1_idx
        error('get_rts: top1_idx inválido (=%g) fuera de rango [1, %d].', top1_idx, numel(nt));
    end

    % 5.5) Obtener cur_pos_top1 = nt(top1_idx) y verificar antes de indexar c_isorts
    cur_pos_top1 = nt(top1_idx);
    if cur_pos_top1 < 1 || cur_pos_top1 > numel(c_isorts) || floor(cur_pos_top1) ~= cur_pos_top1
        error('get_rts: cur_pos_top1 inválido (=%g) fuera de rango [1, %d].', ...
            cur_pos_top1, numel(c_isorts));
    end

    % 5.6) Calcular top1_rt y top1_inten_sum
    top1_rt        = c_isorts(cur_pos_top1);
    top1_inten_sum = inten_sum(top1_idx);

    %% 6) FILTRAR POR top1_rt < 4 minutos
    if top1_rt < 4
        rts            = [];
        top1_rt        = [];
        inten_sum      = [];
        top1_inten_sum = [];
        return;
    end

end  % function get_rts
