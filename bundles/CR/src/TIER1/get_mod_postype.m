function [modpos,modtype] = get_mod_postype(c_sequence,c_modification,Mods)
%%

pos1 = strfind(c_modification,',');
modpos = zeros([1,length(pos1)]);
modtype = zeros([1,length(pos1)]);
if 1==isempty(pos1);
    return;
end;

names = Mods.name;

peplen = length(c_sequence);
pos2 = [0 strfind(c_modification,';')];
for jno = 1 : length(pos1)
    pos = str2num( c_modification(pos2(jno)+1:pos1(jno)-1) );%#ok
    if peplen+1==pos
        pos = peplen;
    end;
    modpos(jno) = pos;
    cmod = c_modification(pos1(jno)+1:pos2(jno+1)-1);
    [tf,loc] = ismember(cmod,names);
    if 1==tf
        modtype(jno) = loc;
    else
        fprintf(1,'no %s in Mods\n',cmod);
        return;
    end;
end;