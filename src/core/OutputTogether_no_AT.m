function OutputTogether_Robust(layout_path, raw_names)
%% OutputTogether_Robust — Consolida muestras con sinónimos y NaN cuando falte
% Corrige mezcla string/char. Usa solo celdas de char. Incluye DEBUG.
%
% Uso:
%   OutputTogether_Robust('E:\EpiProfile_20_AT\src\RawData\histone_layouts', raw_names);

fprintf('[OTR] layout_path = %s\n', layout_path);

%% 0) Sinónimos (canónico -> alternativas en orden de prioridad)
syn.canon   = {'H3_02_9_17','H3_04v3_27_40'};
syn.alts{1} = {'H3_02_9_17','H3_02a_9_17','H3_02b_9_17'};
syn.alts{2} = {'H3_04v3_27_40','H3_04a_27_40'};

    function c = canon_name(basename)
        % basename: char
        c = basename;
        for ii = 1:numel(syn.canon)
            if any(strcmp(basename, syn.alts{ii}))
                c = syn.canon{ii};
                return;
            end
        end
    end

%% 1) Unión de paneles (.mat) canonizados
all_names = {};                           % celda de char sin .mat
per_sample_lists = cell(numel(raw_names),1);

for i = 1:numel(raw_names)
    prefix = sprintf('%02d_%s', i, raw_names{i});
    cur_outpath = fullfile(layout_path, prefix, 'detail');
    D = dir(fullfile(cur_outpath,'*.mat'));
    here = cell(numel(D),1);
    for j=1:numel(D)
        [~, base, ~] = fileparts(D(j).name);  % char
        here{j} = base;
    end
    per_sample_lists{i} = here;
    all_names = [all_names; here]; %#ok<AGROW>
end

% canonizar (celda→celda) y unificar
canon_all = cellfun(@(s) canon_name(s), all_names, 'UniformOutput', false);
canon_set = unique(canon_all);                    % celda de char
canon_set = sort(canon_set);
fprintf('[OTR] Paneles canónicos (%d): %s\n', numel(canon_set), strjoin(canon_set, ', '));

%% 2) Pre-conteo de filas por panel (busca cualquier muestra que lo tenga)
npeps_total = 0;
panel_sizes = containers.Map('KeyType','char','ValueType','double');   % canon -> nrows

for p = 1:numel(canon_set)
    cname = canon_set{p};
    found = false;
    % alternativas de este canon
    alts = {};
    for ii=1:numel(syn.canon)
        if strcmp(cname, syn.canon{ii}), alts = syn.alts{ii}; break; end
    end
    if isempty(alts), alts = {cname}; end

    % busca primera muestra que tenga alguna alternativa
    for i = 1:numel(raw_names)
        prefix = sprintf('%02d_%s', i, raw_names{i});
        cur_outpath = fullfile(layout_path, prefix, 'detail');
        chosen = '';
        for a = 1:numel(alts)
            fp = fullfile(cur_outpath, [alts{a}, '.mat']);
            if exist(fp,'file'), chosen = fp; break; end
        end
        if ~isempty(chosen)
            S = load(chosen, 'His');
            if ~isfield(S,'His')
                warning('[OTR] %s sin His, ignoro: %s', cname, chosen);
                continue;
            end
            nrows = size(S.His.pep_mz,1);
            panel_sizes(cname) = nrows;
            npeps_total = npeps_total + nrows;
            found = true;
            break;
        end
    end
    if ~found
        warning('[OTR] Panel canónico sin datos en ninguna muestra: %s', cname);
        panel_sizes(cname) = 0;
    end
end

%% 3) Matriz grande [npeps_total x (nSamples*3)] (RT, Area, Ratio)
info = nan(npeps_total, numel(raw_names)*3);

%% 4) Volcado por panel canónico y muestra
m = 0;  % desplazamiento global de filas
for p = 1:numel(canon_set)
    cname = canon_set{p};
    nrows = panel_sizes(cname);
    if nrows == 0, continue; end

    fprintf('[OTR] === Panel %s (n=%d) ===\n', cname, nrows);

    % alternativas de este canon
    alts = {};
    for ii=1:numel(syn.canon)
        if strcmp(cname, syn.canon{ii}), alts = syn.alts{ii}; break; end
    end
    if isempty(alts), alts = {cname}; end

    for i = 1:numel(raw_names)
        prefix = sprintf('%02d_%s', i, raw_names{i});
        cur_outpath = fullfile(layout_path, prefix, 'detail');

        % elegir mejor alternativa disponible
        chosen = '';
        chosen_name = '';
        for a = 1:numel(alts)
            fp = fullfile(cur_outpath, [alts{a}, '.mat']);
            if exist(fp,'file'), chosen = fp; chosen_name = alts{a}; break; end
        end

        if isempty(chosen)
            fprintf('[OTR][WARN] %02d:%s → NO hay %s (relleno NaN)\n', i, raw_names{i}, cname);
        else
            S = load(chosen);  % espera 'auc' y 'His'
            if ~isfield(S,'auc')
                fprintf('[OTR][WARN] %02d:%s → %s sin auc (NaN)\n', i, raw_names{i}, chosen_name);
            else
                auc = S.auc;
                block = nan(nrows,3);
                r = min(nrows, size(auc,1));
                c = min(3, size(auc,2));
                block(1:r,1:c) = auc(1:r,1:c);

                % Cols: [(i-1)*3+1 : (i-1)*3+3] = [RT Area Ratio]
                info(m+1:m+nrows, (i-1)*3+1:(i-1)*3+3) = block;
                fprintf('[OTR][OK]  %02d:%s → uso %-15s como %-15s\n', ...
                        i, raw_names{i}, chosen_name, cname);
            end
        end
    end

    m = m + nrows;
end

%% 5) Escribir TSV
[raw_path,fname] = fileparts(layout_path);
if strcmp(fname,'histone_layouts')
    out_tsv = fullfile(raw_path, 'histone_ratios_ROBUST.tsv');
else
    out_tsv = fullfile(fileparts(raw_path), ['histone_ratios_ROBUST_', fname, '.tsv']);
end
fp = fopen(out_tsv,'w');
if fp==-1, error('No puedo abrir salida: %s', out_tsv); end

% Cabeceras (Ratio | Area | RT)
for i=1:numel(raw_names), fprintf(fp,'\t%d,%s', i, raw_names{i}); end; fprintf(fp,'\t');
for i=1:numel(raw_names), fprintf(fp,'\t%d,%s', i, raw_names{i}); end; fprintf(fp,'\t');
for i=1:numel(raw_names), fprintf(fp,'\t%d,%s', i, raw_names{i}); end; fprintf(fp,'\r\n');
fprintf(fp,'Peptide');
for i=1:numel(raw_names), fprintf(fp,'\tRatio'); end; fprintf(fp,'\t');
for i=1:numel(raw_names), fprintf(fp,'\tArea');  end; fprintf(fp,'\t');
for i=1:numel(raw_names), fprintf(fp,'\tRT(min)'); end; fprintf(fp,'\r\n');

% Imprime por grupos usando solo el nombre canónico (estable)
m = 0;
for p = 1:numel(canon_set)
    cname = canon_set{p};
    nrows = panel_sizes(cname);
    if nrows == 0, continue; end
    fprintf(fp, '%s\r\n', cname);
    for r = 1:nrows
        fprintf(fp, '%s_row%02d', cname, r);
        % Bloque Ratio (col 3)
        for i=1:numel(raw_names), fprintf(fp,'\t%f', info(m+r, (i-1)*3+3)); end; fprintf(fp,'\t');
        % Bloque Area (col 2)
        for i=1:numel(raw_names), fprintf(fp,'\t%e', info(m+r, (i-1)*3+2)); end; fprintf(fp,'\t');
        % Bloque RT (col 1)
        for i=1:numel(raw_names), fprintf(fp,'\t%.2f', info(m+r, (i-1)*3+1)); end; fprintf(fp,'\r\n');
    end
    m = m + nrows;
end
fclose(fp);
fprintf('[OTR] Escrito: %s\n', out_tsv);
end
