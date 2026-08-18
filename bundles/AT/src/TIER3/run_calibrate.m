function run_calibrate()
% Script minimo: ejecuta solo DrawISOProfile0 para generar histone_layouts/0_ref_info.mat
% sin correr DrawISOProfile1 (que procesa los 119 modulos por cada raw).
% Util para validar el init_histone0 reparado.

try
    [bOK,raw_path,norganism,nsource,nsubtype] = ReadInput('paras.txt');
    if 0==bOK
        fprintf(1,'FAIL: ReadInput\n');
        return;
    end;
    [def_ptol,soutput,nfigure,ndebug,raw_names] = check_otherparas(raw_path);
    if 1==isempty(raw_names)
        fprintf(1,'FAIL: no raws\n');
        return;
    end;

    % Build special struct (same as EpiProfile.m)
    special.raw_path = raw_path;
    special.nsource = nsource;
    special.nsubtype = nsubtype;
    special.norganism = norganism;
    special.soutput = soutput;
    special.nfigure = nfigure;
    special.ndebug = ndebug;
    if 4==nsource && (0~=nsubtype && 2~=nsubtype)
        special.nhmass = 1;
    else
        special.nhmass = 0;
    end;
    ptol = def_ptol;

    fprintf(1,'raw_path: %s\n',raw_path);
    fprintf(1,'n_raws: %d  ptol: %d\n',length(raw_names),ptol);
    fprintf(1,'calling DrawISOProfile0...\n');
    t1 = clock;
    DrawISOProfile0(raw_path,raw_names,ptol,special);
    t2 = clock;
    fprintf(1,'elapsed: %.1f sec\n',etime(t2,t1));

    mat_file = fullfile(raw_path,'histone_layouts','0_ref_info.mat');
    if 2==exist(mat_file,'file')
        fprintf(1,'OK: %s\n',mat_file);
    else
        fprintf(1,'FAIL: 0_ref_info.mat not created\n');
    end;
catch ME
    fprintf(1,'ERROR: %s\n',ME.message);
    for k=1:length(ME.stack)
        fprintf(1,'  at %s:%d\n',ME.stack(k).name,ME.stack(k).line);
    end;
end;
