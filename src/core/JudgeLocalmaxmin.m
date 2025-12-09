function [nt,nb] = JudgeLocalmaxmin(w)
%% JudgeLocalmaxmin
% Function: get the local max and min indices in a 1D signal/trace.
% INPUT:
%   w : numeric vector of intensities (row or column)
% OUTPUT:
%   nt : indices of local maxima ("tops") in the ORIGINAL indexing of w
%   nb : indices of local minima ("bottoms") in the ORIGINAL indexing of w
%
% Design notes:
% - Edge-aware: treat the first/last points as min/max depending on the
%   initial/final slope.
% - Plateau handling: collapse consecutive equal values before detecting
%   slope changes; then map detected extrema back to original indices.
%
% Revision history:
%   01/19/2010 - the first version (original author)
%   (Comments extensively expanded for clarity; logic unchanged.)

% ---------------------------
% Quick handling of short w's
% ---------------------------
len = length(w);
if len==1
    % Single point: both min and max are trivially at index 1
    nb = 1;
    nt = 1;
    return;
elseif len==2
    % Two points: by definition, the smaller index is the minimum and the
    % larger is the maximum; if equal, set nb=1, nt=2 (arbitrary but stable)
    if w(1)==w(2)
        nb = 1;
        nt = 2;
    elseif w(1)<w(2)
        nb = 1;
        nt = 2;
    else
        nb = 2;
        nt = 1;
    end
    return;
end

% ---------------------------------------------------
% Collapse same-valued adjacents to avoid flat plateaus
% ---------------------------------------------------
% diff(w) == 0 marks positions where w(k) == w(k+1)
x = diff(w);
i = find(x==0);            % indices in 1..len-1 where a plateau occurs
ii = setdiff(1:len,i);     % keep original indices not equal to those 'i'
                           % (this effectively drops one side of each plateau)

% -------------------------------
% Work on the plateau-free vector
% -------------------------------
nw = w(ii);                % filtered series without consecutive duplicates
nlen = length(nw);
if nlen==1
    % Everything collapsed to a single representative value:
    % treat position 1 as min and position len as max in the ORIGINAL space.
    nb = 1;
    nt = len;
    return;
elseif nlen==2
    % After collapsing, just two samples remain; map their roles to the
    % ORIGINAL endpoints for consistency.
    if nw(1)==nw(2)
        nb = 1;
        nt = len;
    elseif nw(1)<nw(2)
        nb = 1;
        nt = len;
    else
        nb = len;
        nt = 1;
    end
    return;
end

% ---------------------------------------------
% Core: detect slope sign changes on plateau-free
% ---------------------------------------------
nx = diff(nw);             % discrete slope between consecutive points
xl = nx(1:end-1);          % left neighbors of the interior points
xr = nx(2:end);            % right neighbors of the interior points
neartop = (xl>0 & xr<0);   % + then -  => local maximum
nearbot = (xl<0 & xr>0);   % - then +  => local minimum
xx = 1:length(nx);         % helper index for nx positions

% ---------------------------
% Build tentative minima list
% ---------------------------
nb = [];
if nx(1)>0
    % If the very first slope is positive, the first point is a local min
    nb(1) = 1;
end
% Interior minima sit at transitions - to +
nb = [nb xx(nearbot)+1];   % +1 to move from slope index to point index
if nx(end)<0
    % If the last slope is negative, the last point is a local min
    nb = [nb nlen];
end

% ---------------------------
% Build tentative maxima list
% ---------------------------
nt = [];
if nx(1)<0
    % If the very first slope is negative, the first point is a local max
    nt(1) = 1;
end
% Interior maxima sit at transitions + to -
nt = [nt xx(neartop)+1];   % +1 to move from slope index to point index
if nx(end)>0
    % If the last slope is positive, the last point is a local max
    nt = [nt nlen];
end

% ------------------------------------
% Map extrema back to ORIGINAL indices
% ------------------------------------
% We detected extrema on 'nw' (plateau-free). 'ii' holds the kept original
% indices, so use it to translate the positions back to the original w.
nb = ii(nb);
nt = ii(nt);

% ---------------------------------------------------------
% Final edge consistency check in ORIGINAL index coordinates
% ---------------------------------------------------------
% Ensure the first extrema align with boundaries when appropriate.
if nt(1)<nb(1) && nt(1)>1
    % If the first max is before the first min and it's not at index 1,
    % force the first max to be the boundary (index 1).
    nt(1) = 1;
elseif nb(1)<nt(1) && nb(1)>1
    % Symmetrically for the first min
    nb(1) = 1;
end

% Ensure the last extrema align with the end boundary when appropriate.
if nt(end)>nb(end) && nt(end)<len
    % Last max comes after last min but not at the very end: pin to 'len'.
    nt(end) = len;
elseif nb(end)>nt(end) && nb(end)<len
    % Symmetrically for the last min
    nb(end) = len;
end
