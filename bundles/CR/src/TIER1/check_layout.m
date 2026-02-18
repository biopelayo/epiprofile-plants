function check_layout()
% the cursor goes to the last line, and press F12 to set a breakpoint
% press F5 to start, and check layouts for RAW
% if saving changes, press F5; otherwise, press red button above 'Quit Debugging'

clc;

[bOK,raw_path] = ReadInput('paras.txt');% read paras
if 0==bOK
    return;
end;

mat_file = fullfile(raw_path,'histone_layouts','0_ref_info.mat');
if 0==exist(mat_file,'file')
    return;
end;
load(mat_file);% load ref

fprintf(1,'checking layouts...\n');
save(mat_file,'AllUnHis');% save changes