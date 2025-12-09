function [cur_rts, cur_intens, cur_mono_isointens] = get_histone11(MS1_index, MS1_peaks, ptol, unitdiff, His, hno)
%%
% GET_HISTONE11 Extrae los valores de retención (RT) y las intensidades monoisotópicas
% de un péptido/histona dado su número de carga. Trabaja sobre un intervalo de retención
% alrededor de His.rt_ref(hno) ± 0.5 minutos.  
%
%   [cur_rts, cur_intens, cur_mono_isointens] = get_histone11(MS1_index, MS1_peaks, ptol, unitdiff, His, hno)
%
%   ENTRADAS:
%     • MS1_index   : Matriz N×M donde la segunda columna (MS1_index(:,2)) contiene
%                     los tiempos de retención (RT) de cada escaneo MS1.
%     • MS1_peaks   : Matriz que contiene las intensidades asociadas a cada línea/scan
%                     de MS1_index.
%     • ptol        : Tolerancia en m/z (en unidades de ppm o unidades propias) para 
%                     la extracción de perfil de isótopos. Si se provee como 100, internamente
%                     se ajusta a 10.
%     • unitdiff    : Separación en masa isotópica (aproximadamente 1 Da, pero ajustada
%                     por la carga).
%     • His         : Estructura con los metadatos de la histona/péptido:
%         – His.pep_mz    : Matriz [npep × ncharge] con valores de m/z teóricos de cada
%                           péptido/histona para cada número de carga.
%         – His.pep_ch    : Matriz [npep × ncharge] con cada número de carga asociado.
%         – His.rt_ref    : Vector [npep × 1] con el tiempo de retención de referencia de
%                           cada péptido. Si algún elemento es 0 o negativo, se considera
%                           “no estimado” y la función devolverá ceros.
%     • hno         : Índice (entero) del péptido/histona que se está procesando.
%
%   SALIDAS:
%     • cur_rts            : Vector [1 × ncharge] que contiene la RT encontrada (posiblemente
%                            ajustada al máximo pico cercano) para cada número de carga.
%     • cur_intens         : Vector [1 × ncharge] con la intensidad (área integrada) para
%                            cada número de carga. Si no se detecta pico en ±0.5 min, se deja 0.
%     • cur_mono_isointens : Vector columna con la intensidad monoisotópica (isótopo j=2)
%                            extraída para la carga 1 (solo se conserva en la primera iteración).
%
%   NOTA SOBRE EL ERROR “Array indices must be positive integers or logical values”:
%     Este error ocurre cuando MATLAB intenta indexar un vector o matriz usando un índice que
%     no es un entero positivo válido (por ejemplo, 0, negativo, decimal, vacío o mayor que el
%     tamaño del arreglo). En este código, varias líneas hacen indexing directo (p(1), pp(end),
%     c_isorts(cur_pos), nb(x(end)+1), etc.) y todas esas operaciones asumen que las variables
%     involucradas no están vacías y contienen enteros válidos.  
%
%     A continuación se documenta cada paso y se introducen mejoras (validaciones) para:
%       1) No procesar hno cuyo His.rt_ref(hno) = 0, porque la ventana [rt_ref-0.5, rt_ref+0.5]
%          no intercepta ningún escaneo MS1 real (RTs típicos comienzan en ~1).  
%       2) Comprobar que los vectores p y pp (de índices de MS1_index(:,2)) no estén vacíos antes 
%          de extraer p(1) o pp(end).  
%       3) Validar que cur_pos (extraído de nt) sea un entero en rango antes de indexar c_isorts.  
%       4) Validar que nb(x(end)+1) esté dentro de límites antes de indexar nb.  
%
%   USO TÍPICO:
%     valid_hno = find(His.rt_ref > 0);
%     for k = 1:numel(valid_hno)
%         hno = valid_hno(k);
%         [rts, intens, mono] = get_histone11(MS1_index, MS1_peaks, ptol, unitdiff, His, hno);
%         % … almacenar resultados …
%     end
%
%   Si se llama directamente con un hno tal que His.rt_ref(hno) == 0, la función emite un warning
%   y devuelve vectores de ceros, evitando así errores de índices vacíos.
%

    %% 1) Verificación inicial de que rt_ref(hno) esté estimado (> 0)
    if His.rt_ref(hno) <= 0
        warning( ...
            'get_histone11: rt_ref(%d) = %.4g. No se extrae perfil MS1 (se devuelven ceros).', ...
            hno, His.rt_ref(hno) ...
        );
        % Determinar ncharge para inicializar los vectores de salida
        [~, ncharge] = size(His.pep_mz);  %#ok
        cur_rts            = zeros(1, ncharge);
        cur_intens         = zeros(1, ncharge);
        cur_mono_isointens = [];  % No se extrae ninguna serie monoisotópica
        return
    end

    %% 2) Inicialización de vectores de salida a cero (por defecto)
    [~, ncharge] = size(His.pep_mz);  %#ok
    cur_rts    = zeros([1, ncharge]);
    cur_intens = zeros([1, ncharge]);
    % cur_mono_isointens solo se asignará a la carga 1 en el bucle

    %% 3) Definir la ventana de retención (RT) alrededor de rt_ref(hno)
    delta  = 0.5;  % ±0.5 minutos de tolerancia de RT
    rt_min = His.rt_ref(hno) - delta;
    rt_max = His.rt_ref(hno) + delta;

    % 3.1) Encontrar posición de inicio en MS1_index(:,2) >= rt_min
    p = find(MS1_index(:,2) >= rt_min);
    if isempty(p)
        % Ningún escaneo MS1 con RT >= rt_min -> no se puede extraer perfil
        error( ...
            'get_histone11: No hay ningún RT_MS1 >= %.4g (rt_ref - %.4g).', ...
            rt_min, delta ...
        );
    end
    rt_i1 = p(1);  % Primer índice válido (el más cercano a rt_min)

    % 3.2) Encontrar posición de fin en MS1_index(:,2) <= rt_max
    pp = find(MS1_index(:,2) <= rt_max);
    if isempty(pp)
        % Ningún escaneo MS1 con RT <= rt_max -> no se puede extraer perfil
        error( ...
            'get_histone11: No hay ningún RT_MS1 <= %.4g (rt_ref + %.4g).', ...
            rt_max, delta ...
        );
    end
    rt_i2 = pp(end);  % Último índice válido (el más cercano a rt_max)

    %% 4) Ajuste especial de ptol
    % En algunos casos ptol se pasa como 100 para marcar un valor “por defecto”.
    % Internamente, si ptol==100, lo cambiamos a 10 (valor real de tolerancia).
    if ptol == 100
        ptol = 10;
    end

    %% 5) Bucle sobre cada número de carga (1 a ncharge)
    for jno = 1:ncharge
        %------------------------------------------------------------
        % 5.1) Extraer m/z teóricos e isótopos separados por unidad de masa
        %------------------------------------------------------------
        c_mz = His.pep_mz(hno, jno);  % m/z monoisotópico teórico
        c_ch = His.pep_ch(hno, jno);  % número de carga

        % Construir vector de m/z de isótopos que vamos a extraer:
        %   1) monoisótopo:   c_mz
        %   2) isótopo -1:    c_mz - (unitdiff / c_ch)
        %   3) isótopo +1:    c_mz + (unitdiff / c_ch)
        %   4) isótopo +2:    c_mz + 2*(unitdiff / c_ch)
        c_ref_isomzs = [ ...
            c_mz - unitdiff / c_ch, ...
            c_mz, ...
            c_mz + unitdiff / c_ch, ...
            c_mz + 2 * unitdiff / c_ch ...
        ];

        % Determinar si consideramos el isótopo +1 (nC13 = 1) o no (nC13 = 0)
        if ptol > 100 && c_ch >= 3
            nC13 = 1;
        else
            nC13 = 0;
        end

        %------------------------------------------------------------
        % 5.2) Llamar a GetProfiles para extraer el perfil MS1 entre rt_i1 y rt_i2
        %------------------------------------------------------------
        %    [c_isorts, c_ref_isointens] = GetProfiles(...)
        %
        %    c_isorts        : Vector columna con tiempos de retención (ordenados)
        %                       correspondientes a escaneos en el rango [rt_i1:rt_i2].
        %    c_ref_isointens : Matriz [num_steps × 4] con las intensidades en cada RT
        %                       para cada uno de los 4 m/z de isótopos.
        [c_isorts, c_ref_isointens] = GetProfiles( ...
            MS1_index, MS1_peaks, c_ref_isomzs, c_ch, ptol, nC13, rt_i1:rt_i2 ...
        );

        % Tomar la columna j=2 de c_ref_isointens, que corresponde al monoisótopo
        j = 2;
        c_mono_isointens = c_ref_isointens(:, j);

        % Guardar el vector de intensidades monoisotópicas completo para la carga 1
        if jno == 1
            cur_mono_isointens = c_mono_isointens;
        end

        %------------------------------------------------------------
        % 5.3) Llamar a GetTopBottom11 para identificar “tops” y “bottoms”
        %------------------------------------------------------------
        %    [nt, nb, top1_idx, inten_sum] = GetTopBottom11(c_mono_isointens)
        %
        %    nt       : Vector con índices (relativos a c_mono_isointens) donde hay picos “top”.
        %    nb       : Vector con índices (relativos a c_mono_isointens) donde hay picos “bottom”.
        %    top1_idx : Órbita (índice dentro de nt) del pico candidato con mayor área parcial.
        %    inten_sum: Suma de intensidades en torno a cada posición de nt.
        [nt, nb, top1_idx, inten_sum] = GetTopBottom11(c_mono_isointens); %#ok

        %------------------------------------------------------------
        % 5.4) Filtrar solo aquellos “tops” cuya RT (c_isorts(nt(i))) esté a ±0.5 de rt_ref
        %------------------------------------------------------------
        flag = zeros([1, length(nt)]);  % Inicializo máscara de ceros
        for i = 1:length(nt)
            % nt(i) es un índice en c_isorts; verifico si su RT está dentro de la ventana
            if abs(c_isorts(nt(i)) - His.rt_ref(hno)) <= 0.5
                flag(i) = 1;
            end
        end
        x = find(flag == 1);  % Índices de nt que cumplen la condición

        %------------------------------------------------------------
        % 5.5) Si existe al menos un “top” cercano a rt_ref, calcular RT e intensidad
        %------------------------------------------------------------
        if ~isempty(x)
            % 5.5.1) Elegir el pico con mayor inten_sum entre los candidatos
            [~, id] = max(inten_sum(x));  %#ok
            top1_idx = x(id);             % Posición en nt del pico elegido

            % 5.5.2) Obtener el índice real en c_isorts
            cur_pos = nt(top1_idx);

            % 5.5.3) Validar que cur_pos sea un entero válido dentro de [1, numel(c_isorts)]
            if floor(cur_pos) ~= cur_pos || cur_pos < 1 || cur_pos > numel(c_isorts)
                % Si no es un índice válido, asigno valor por defecto (RT_ref e intensidad 0)
                warning( ...
                    'get_histone11: cur_pos inválido (=%g) en hno=%d, jno=%d. Intensidad=0.', ...
                    cur_pos, hno, jno ...
                );
                cur_rts(jno)    = His.rt_ref(hno);
                cur_intens(jno) = 0;
            else
                % 5.5.4) Si es válido, extraigo la RT real del pico
                cur_rts(jno) = c_isorts(cur_pos);

                % 5.5.5) Si hay varios candidatos en x, ajusto nb para el cálculo de get_area
                if length(x) >= 2
                    idx1 = x(1);
                    idx2 = x(end) + 1;
                    % Verificar que idx2 no exceda la longitud de nb
                    if idx2 <= length(nb)
                        nb = [nb(idx1), nb(idx2)];
                    else
                        warning( ...
                            'get_histone11: nb(x(end)+1) fuera de rango para hno=%d, jno=%d.', ...
                            hno, jno ...
                        );
                        % Mantener nb original (podría usarse un valor por defecto)
                    end
                end

                % 5.5.6) Calcular el área (intensidad integrada) alrededor de cur_pos
                cur_intens(jno) = get_area( ...
                    c_isorts, c_ref_isointens, nb, cur_pos, c_mz, c_ch, ...
                    MS1_index, MS1_peaks, unitdiff, ptol ...
                );
            end

        else
            %------------------------------------------------------------
            % 5.6) Si no hay ningún pico dentro de ±0.5 de rt_ref, asigno RT_ref y 0
            %------------------------------------------------------------
            cur_rts(jno)    = His.rt_ref(hno);
            cur_intens(jno) = 0;
        end

    end  % for jno = 1:ncharge

end  % function get_histone11
