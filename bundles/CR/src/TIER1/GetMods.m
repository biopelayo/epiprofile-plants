function Mods = GetMods()
%% Function: get the mod table
% Output:
%   Mods - the mod table
% http://www.unimod.org/

bmono = 1;
bavge = 2;
nMod = 500;% initial number of modifications
n = 0;% real number of modifications

% premalloc for mass
Mods.mass = zeros([2,nMod]);

n = n + 1;
Mods.set{n} = 'KSTY';
Mods.name{n} = 'pr';
Mods.mass(bmono,n) = 56.026215;
Mods.mass(bavge,n)  = 56.0633;

n = n + 1;
Mods.set{n} = 'KSTY';
Mods.name{n} = 'ac';
Mods.mass(bmono,n) = 42.010565;
Mods.mass(bavge,n)  = 42.0367;

n = n + 1;
Mods.set{n} = 'K';
Mods.name{n} = 'me1';
Mods.mass(bmono,n) = 70.041865;
Mods.mass(bavge,n)  = 70.0898;

n = n + 1;
Mods.set{n} = 'R';
Mods.name{n} = 'Me1';
Mods.mass(bmono,n) = 14.015650;
Mods.mass(bavge,n)  = 14.0266;

n = n + 1;
Mods.set{n} = 'KR';
Mods.name{n} = 'me2';
Mods.mass(bmono,n) = 28.031300;
Mods.mass(bavge,n)  = 28.0532;

n = n + 1;
Mods.set{n} = 'K';
Mods.name{n} = 'me3';
Mods.mass(bmono,n) = 42.046950;
Mods.mass(bavge,n)  = 42.0797;

n = n + 1;
Mods.set{n} = 'STY';
Mods.name{n} = 'ph';
Mods.mass(bmono,n) = 79.966331;
Mods.mass(bavge,n)  = 79.9799;

n = n + 1;
Mods.set{n} = 'K';
Mods.name{n} = 'cr';
Mods.mass(bmono,n) = 68.026215;
Mods.mass(bavge,n)  = 68.0740;

n = n + 1;
Mods.set{n} = 'K';
Mods.name{n} = 'ub';
Mods.mass(bmono,n) = 114.042927;
Mods.mass(bavge,n)  = 114.1026;

n = n + 1;
Mods.set{n} = 'M';
Mods.name{n} = 'ox';
Mods.mass(bmono,n) = 15.994915;
Mods.mass(bavge,n)  = 15.9994;

n = n + 1;
Mods.set{n} = 'K';
Mods.name{n} = 'lak';
Mods.mass(bmono,n) = 8.014199;
Mods.mass(bavge,n)  = 7.9427;

n = n + 1;
Mods.set{n} = 'R';
Mods.name{n} = 'lar';
Mods.mass(bmono,n) = 10.008269;
Mods.mass(bavge,n)  = 9.9296;

n = n + 1;
Mods.set{n} = 'K';
Mods.name{n} = 'hac';
Mods.mass(bmono,n) = 44.017274;
Mods.mass(bavge,n)  = 44.0220;

n = n + 1;
Mods.set{n} = 'K';
Mods.name{n} = 'hme1';
Mods.mass(bmono,n) = 73.060695;
Mods.mass(bavge,n)  = 73.1083;

n = n + 1;
Mods.set{n} = 'K';
Mods.name{n} = 'hme2';
Mods.mass(bmono,n) = 34.06896;
Mods.mass(bavge,n)  = 34.0902;

n = n + 1;
Mods.set{n} = 'K';
Mods.name{n} = 'hme3';
Mods.mass(bmono,n) = 51.10344;
Mods.mass(bavge,n)  = 51.1352;

n = n + 1;
Mods.set{n} = 'K';
Mods.name{n} = 'prlak';
Mods.mass(bmono,n) = 64.040414;
Mods.mass(bavge,n)  = 64.006;

n = n + 1;
Mods.set{n} = 'K';
Mods.name{n} = 'aclak';
Mods.mass(bmono,n) = 50.024764;
Mods.mass(bavge,n)  = 49.9794;

n = n + 1;
Mods.set{n} = 'K';
Mods.name{n} = 'me1lak';
Mods.mass(bmono,n) = 78.056064;
Mods.mass(bavge,n)  = 78.0325;

n = n + 1;
Mods.set{n} = 'K';
Mods.name{n} = 'me2lak';
Mods.mass(bmono,n) = 36.045499;
Mods.mass(bavge,n)  = 35.9959;

n = n + 1;
Mods.set{n} = 'K';
Mods.name{n} = 'me3lak';
Mods.mass(bmono,n) = 50.061149;
Mods.mass(bavge,n)  = 50.0224;

n = n + 1;
Mods.set{n} = 'K';
Mods.name{n} = 'me11';
Mods.mass(bmono,n) = 74.06405;
Mods.mass(bavge,n)  = 74.1009;

n = n + 1;
Mods.set{n} = 'K';
Mods.name{n} = 'me21';
Mods.mass(bmono,n) = 32.053485;
Mods.mass(bavge,n)  = 32.0643;

n = n + 1;
Mods.set{n} = 'K';
Mods.name{n} = 'me22';
Mods.mass(bmono,n) = 36.07567;
Mods.mass(bavge,n)  = 36.0754;

n = n + 1;
Mods.set{n} = 'K';
Mods.name{n} = 'me31';
Mods.mass(bmono,n) = 46.069135;
Mods.mass(bavge,n)  = 46.0908;

n = n + 1;
Mods.set{n} = 'K';
Mods.name{n} = 'me32';
Mods.mass(bmono,n) = 50.09132;
Mods.mass(bavge,n)  = 50.1019;

n = n + 1;
Mods.set{n} = 'K';
Mods.name{n} = 'me33';
Mods.mass(bmono,n) = 54.113505;
Mods.mass(bavge,n)  = 54.113;

n = n + 1;
Mods.set{n} = 'M';
Mods.name{n} = 'cd3';
Mods.mass(bmono,n) = 4.022185;
Mods.mass(bavge,n)  = 4.0111;

if n<nMod
    Mods.mass = Mods.mass(1:2,1:n);
end;