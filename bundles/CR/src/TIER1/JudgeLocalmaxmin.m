function [nt,nb] = JudgeLocalmaxmin(w)
%% JudgeLocalmaxmin
% Function: get the local max and min
% Inputs:
%    w - the intensity values
% Outputs:
%    nt - the localmax
%    nb - the localmin
%
% Revision history:
%    01/19/2010 - the first version

% check the number of w
len = length(w);
if len==1
    nb = 1;
    nt = 1;
    return;
elseif len==2
    if w(1)==w(2)
        nb = 1;
        nt = 2;
    elseif w(1)<w(2)
        nb = 1;
        nt = 2;
    else
        nb = 2;
        nt = 1;
    end;
    return;
end;

% cut off the same adjacent values
x = diff(w);
i = find(x==0);
ii = setdiff(1:len,i);

% check the number of nw 
nw = w(ii);
nlen = length(nw);
if nlen==1
    nb = 1;
    nt = len;
    return;
elseif nlen==2
    if nw(1)==nw(2)
        nb = 1;
        nt = len;
    elseif nw(1)<nw(2)
        nb = 1;
        nt = len;
    else
        nb = len;
        nt = 1;
    end;
    return;
end;

% get nb and nt
nx = diff(nw);
xl = nx(1:end-1);% left side
xr = nx(2:end);	  % right side
neartop = (xl>0 & xr<0);
nearbot = (xl<0 & xr>0);
xx = 1:length(nx);

nb = [];
if nx(1)>0
    nb(1) = 1;
end;
nb = [nb xx(nearbot)+1];% tentative minimum
if nx(end)<0
    nb = [nb nlen];
end;

nt = [];
if nx(1)<0
    nt(1) = 1;
end;
nt = [nt xx(neartop)+1];% tentative maximum
if nx(end)>0
    nt = [nt nlen];
end;

% return to the index of w
nb = ii(nb);
nt = ii(nt);

% check the end of nb and nt
if nt(1)<nb(1) && nt(1)>1
    nt(1) = 1;
elseif nb(1)<nt(1) && nb(1)>1
    nb(1) = 1;
end;

if nt(end)>nb(end) && nt(end)<len
    nt(end) = len;
elseif nb(end)>nt(end) && nb(end)<len
    nb(end) = len;
end;