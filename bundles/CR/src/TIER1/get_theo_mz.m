function [theo_mz,theo_tp] = get_theo_mz(c_sequence,c_modification,chg,Mods,ActiveType)
%%

% init
element = [12 1.0078246 14.0030732 15.9949141 31.972070];% element mass
mH = element(2);
mCO = element(1) + element(4);%#ok
mHO = element(2) + element(4);
mH2O = element(2)*2 + element(4);
mNH3 = element(2)*3 + element(3);
mNH2 = element(2)*2 + element(3);

% peptide
aamass = Getaamass();
idx = c_sequence-'A'+1;
residuemass = aamass(idx,1)';

% modification
deltam = get_mod_mass(c_sequence,c_modification,Mods);

% peptide+modification
residuemass = residuemass + deltam;

% ions and charge
pmass = 1.007276;
if 1==strcmp(ActiveType,'CID')
    Mb = cumsum(residuemass(1:end-1));
    My = cumsum(rot90(rot90(residuemass(2:end))))+mH2O;

    M = [Mb My];
    TP = [ones(1,length(Mb))*1000+(1:length(Mb))*10 2*ones(1,length(My))*1000+(1:length(My))*10];
    [M,ix] = unique(M);
    TP = TP(ix);

    theo_mz = repmat(0,[1,length(M)*chg]);
    theo_tp = repmat(0,[1,length(M)*chg]);
    for ino = 1 : chg
        theo_mz((ino-1)*length(M)+1:ino*length(M)) = (M+ino*pmass)/ino;
        theo_tp((ino-1)*length(M)+1:ino*length(M)) = TP + ino;
    end;
    [theo_mz,ix] = unique(theo_mz);
    theo_tp = theo_tp(ix);
elseif 1==strcmp(ActiveType,'ETD')
    Mc = cumsum(residuemass(1:end-1))+mNH3;
    Mz = cumsum(rot90(rot90(residuemass(2:end))))+mHO-mNH2+mH;

    M = [Mc Mz];
    TP = [ones(1,length(Mc))*1000+(1:length(Mc))*10 2*ones(1,length(Mz))*1000+(1:length(Mz))*10];
    [M,ix] = unique(M);
    TP = TP(ix);

    theo_mz = repmat(0,[1,length(M)*chg]);
    theo_tp = repmat(0,[1,length(M)*chg]);
    for ino = 1 : chg
        theo_mz((ino-1)*length(M)+1:ino*length(M)) = (M+ino*pmass)/ino;
        theo_tp((ino-1)*length(M)+1:ino*length(M)) = TP + ino;
    end;
    [theo_mz,ix] = unique(theo_mz);
    theo_tp = theo_tp(ix);
else
    theo_mz = [];
    theo_tp =  [];
end;