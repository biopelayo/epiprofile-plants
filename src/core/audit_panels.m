%% === AUDITORÍA AUTOCONTENIDA: plantilla vs cada muestra ===
% 1) Ruta base (ajústala si cambia)
layout_path = 'E:\EpiProfile_20_AT\src\RawData\histone_layouts';

% 2) Construir lista de muestras a partir de carpetas NN_<rawname>
d = dir(layout_path);
names = {d([d.isdir]).name};
tokens = regexp(names, '^(\d{2})_(.+)$', 'tokens', 'once');   % { {NN, raw} | [] | ... }
mask   = ~cellfun('isempty', tokens);
nums   = cellfun(@(t) str2double(t{1}), tokens(mask));
rns    = cellfun(@(t) t{2},          tokens(mask), 'UniformOutput', false);
[nums,ord] = sort(nums);                      % orden por NN
raw_names  = rns(ord);                        % nombres ordenados

fprintf('Detectadas %d muestras: %s\n', numel(raw_names), strjoin(raw_names, ', '));

% 3) Plantilla: lista de .mat en la muestra 1
p1 = fullfile(layout_path, sprintf('%02d_%s', nums(1), raw_names{1}), 'detail');
if ~isfolder(p1), error('No existe carpeta: %s', p1); end
templ = dir(fullfile(p1,'*.mat'));
templ_names = sort({templ.name});
fprintf('Plantilla muestra %02d (%s): %d .mat\n', nums(1), raw_names{1}, numel(templ_names));

need_02a = ismember('H3_02a_9_17.mat', templ_names);  % ¿la plantilla espera 02a?

% 4) Comparar cada muestra contra la plantilla
for k = 1:numel(raw_names)
    pk = fullfile(layout_path, sprintf('%02d_%s', nums(k), raw_names{k}), 'detail');
    here = dir(fullfile(pk,'*.mat'));
    here_names = sort({here.name});
    missing = setdiff(templ_names, here_names);
    extra   = setdiff(here_names, templ_names);

    fprintf('Muestra %02d (%s): faltan=%d extra=%d\n', nums(k), raw_names{k}, numel(missing), numel(extra));
    if ~isempty(missing)
        fprintf('   Faltan (hasta 10): %s\n', strjoin(missing(1:min(end,10)), ', '));
    end
    if ~isempty(extra)
        fprintf('   Extra (hasta 10): %s\n', strjoin(extra(1:min(end,10)), ', '));
    end
    if need_02a && ~ismember('H3_02a_9_17.mat', here_names)
        fprintf('   ALERTA: falta H3_02a_9_17.mat en esta muestra\n');
    end
end
