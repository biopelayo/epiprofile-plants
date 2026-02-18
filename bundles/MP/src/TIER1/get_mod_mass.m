function deltam = get_mod_mass(c_sequence,c_modification,Mods)
%%

peplen = length(c_sequence);
deltam = repmat(0,[1,peplen]);
pos1 = strfind(c_modification,',');
if 1==isempty(pos1);
    return;
end;

masses = Mods.mass(1,:);
names = Mods.name;

pos2 = [0 strfind(c_modification,';')];
for jno = 1 : length(pos1)
    pos = str2num( c_modification(pos2(jno)+1:pos1(jno)-1) );%#ok
    if 0==pos
        pos = 1;
    elseif peplen+1==pos
        pos = peplen;
    end;
    cmod = c_modification(pos1(jno)+1:pos2(jno+1)-1);
    [tf,loc] = ismember(cmod,names);
    if 1==tf
        deltam(pos) = deltam(pos) + masses(loc);
    else
        deltam = repmat(0,[1,peplen]);
        fprintf(1,'no %s in Mods\n',cmod);
        return;
    end;
end;