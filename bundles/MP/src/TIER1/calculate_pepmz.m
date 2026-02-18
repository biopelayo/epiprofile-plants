function pep_mz = calculate_pepmz(His)
%%

pep_mz = repmat(0,[size(His.mod_type,1),size(His.pep_ch,2)]);

Mods = GetMods();
aamass = Getaamass();
element = [12 1.0078246 14.0030732 15.9949141 31.972070];% element mass
mH2O = element(2)*2 + element(4);
pmass = 1.007276;

if 0==strcmp(His.pep_seq,'unmod')
    c_seq = His.pep_seq;
    idx = c_seq-'A'+1;
    residuemass = aamass(idx,1)';
    for ino=1:size(His.mod_type,1)
        c_mod = His.mod_type{ino};
        deltam = get_mod_mass(c_seq,c_mod,Mods);
        % peptide+modification
        residuemass_new = residuemass + deltam;
        Mr = sum(residuemass_new)+mH2O;
        for jno=1:size(His.pep_ch,2)
            c_ch = His.pep_ch(ino,jno);
            pep_mz(ino,jno) = (Mr+c_ch*pmass)/c_ch;
        end;
    end;
else
    for ino=1:size(His.mod_type,1)
        c_seq = His.mod_short{ino};
        p = strfind(c_seq,'.');
        if 0==isempty(p)
            c_seq = c_seq(p(end)+1:end);
        end;
        idx = c_seq-'A'+1;
        residuemass = aamass(idx,1)';

        c_mod = His.mod_type{ino};
        deltam = get_mod_mass(c_seq,c_mod,Mods);
        % peptide+modification
        residuemass_new = residuemass + deltam;
        Mr = sum(residuemass_new)+mH2O;
        for jno=1:size(His.pep_ch,2)
            c_ch = His.pep_ch(ino,jno);
            pep_mz(ino,jno) = (Mr+c_ch*pmass)/c_ch;
        end;
    end;
end;