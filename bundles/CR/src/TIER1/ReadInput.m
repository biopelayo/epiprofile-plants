function [bOK,raw_path,norganism,nsource,nsubtype] = ReadInput(para_txtfile)
%%

bOK = 0;
raw_path = '';% the datapath of raw files
norganism = 1;% 1: Human, 2: Mouse
nsource = 1;% 1: histone_normal, 2: histone_SILAC, 3: histone_C13, 4: histone_N15
nsubtype = 0;% if histone_N15, 0: N14 light Mods, 1: N15 light Mods, 2: N14 heavy Mods, 3: N15 heavy Mods, 4: 0+1, 5: 0+3

% check
if 0==exist( para_txtfile,'file' )
    disp(['the para file does not exist: ',para_txtfile]);
    return;
end;

% open
fp = fopen( para_txtfile,'r' );
if -1==fp
    disp(['can not open: ',para_txtfile]);
    return;
end;

% read
while 0==feof(fp)
    str = fgetl(fp);
    % raw_path
    pos = strfind(str,'raw_path');
    while 0==feof(fp) && (1==isempty(pos) || 1~=pos(1))
        str = fgetl(fp);
        pos = strfind(str,'raw_path');
    end;
    if 1==feof(fp)
        break;
    end;
    pos = strfind(str,'=');
    raw_path = str(pos+1:end);

    % norganism
    pos = strfind(str,'norganism');
    while 0==feof(fp) && 1==isempty(pos)
        str = fgetl(fp);
        pos = strfind(str,'norganism');
    end;
    if 1==isempty(pos)
        break;
    end;
    pos = strfind(str,'=');
    norganism = str2double( str(pos+1:end) );

    % nsource
    pos = strfind(str,'nsource');
    while 0==feof(fp) && 1==isempty(pos)
        str = fgetl(fp);
        pos = strfind(str,'nsource');
    end;
    if 1==isempty(pos)
        break;
    end;
    pos = strfind(str,'=');
    nsource = str2double( str(pos+1:end) );

    % nsubtype
    pos = strfind(str,'nsubtype');
    while 0==feof(fp) && 1==isempty(pos)
        str = fgetl(fp);
        pos = strfind(str,'nsubtype');
    end;
    if 1==isempty(pos)
        break;
    end;
    pos = strfind(str,'=');
    nsubtype = str2double( str(pos+1:end) );
end;
fclose(fp);

bOK = 1;