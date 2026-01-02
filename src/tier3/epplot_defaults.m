function style = epplot_defaults(varargin)
% EPPLOT_DEFAULTS — Estilo global coherente para todas las figuras.
% Uso: style = epplot_defaults('FontSize',11,'DPI',600);

    p = inputParser;
    addParameter(p,'FontName','Arial');
    addParameter(p,'FontSize',10);
    addParameter(p,'LineWidth',1.3);
    addParameter(p,'Alpha',0.8);
    addParameter(p,'DPI',600);
    addParameter(p,'FigureSize',[6 4]);      % pulgadas [w h]
    addParameter(p,'WideSize',[10 4.5]);     % pulgadas (overlays)
    addParameter(p,'HeatmapSize',[8 8]);     % pulgadas (heatmaps)
    addParameter(p,'FacetCols',2);           % nº columnas para facet
    parse(p,varargin{:});
    style = p.Results;

    % Paletas seguras y consistentes
    style.Palette.Groups = [ ...
        0.000,0.447,0.741;  % azul
        0.902,0.624,0.000;  % naranja
        0.000,0.620,0.451;  % verdoso
        0.835,0.369,0.000;  % rojo/naranja
        0.800,0.475,0.655;  % malva
        0.337,0.706,0.914;  % azul claro
        0.941,0.894,0.259;  % amarillo
        0.000,0.000,0.000]; % negro
    style.Palette.Divergent = local_diverging_balance(); % para correlaciones/ratios
    style.Palette.Grey = [0.95 0.95 0.95];
end

% ============== helpers locales (quedan dentro del mismo fichero) ==============
function C = local_diverging_balance()
% 256×3 RGB diverging “balance” (−1↔azul, 0↔gris, +1↔rojo), suave y simétrica.
    n = 256; x = linspace(-1,1,n)';
    r = 0.5+0.5*tanh( 2.2*x);
    b = 0.5+0.5*tanh(-2.2*x);
    g = 0.5 + 0.4*(1-abs(x)).^1.3;
    C = [b g r]; C = max(min(C,1),0);
end
