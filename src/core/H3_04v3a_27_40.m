function H3_04v3a_27_40(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,cur_outpath,special)
%%
% H3_04v3a_27_40
% Targeted panel for H3 peptide KSAPSTGGVKKPHR (positions 27–40), i.e. the
% variant where the motif is KSAPST (v3) instead of KSAPAT. This panel
% quantifies S28-bearing states (ph/ac/pr) in combination with the K27
% methylation ladder (me1/me2/me3) and selected K36 states (me1/me2).
%
% Design choices:
% - No anchoring on the unmodified peptide here. Reference RTs are imported
%   from a base panel XLS (H3_04v3_27_40.xls, [rt] section) and used to
%   shape search windows.
% - Two relocation modes:
%     * MS1-only heuristic windows via get_rts()  (special.nDAmode ~= 2)
%     * MS2-assisted windows via get_rts2()       (special.nDAmode == 2)
% - Final per-state extraction uses get_histone12() (consistent S28 extractor).
%
% Inputs:
%   MS1_index, MS1_peaks : MS1 centroids and (scan, RT[min]) index
%   MS2_index, MS2_peaks : MS2 centroids and (scan, RT[min]) index
%   ptol                 : precursor ppm tolerance
%   cur_outpath          : output directory
%   special              : struct with fields
%       raw_path (string) : kept for API compatibility
%       ndebug (0/1)      : debug relocation using relocateD (MS1-only)
%       nDAmode (0/1/2)   : 2 => MS2-assisted relocation; else MS1
%       nhmass            : optional MS2 diagnostic configuration
%
% Outputs (files):
%   <cur_outpath>/H3_04v3a_27_40.mat  : quantitative matrices (pep_rts, pep_intens)
%   Figures via draw_layout()          : XICs, apexes, optional MS2 overlays
%   PSMs via GetPSM() if nDAmode==1    : identification table

% -------------------------------------------------------------------------
% Idempotency check (skip if output already exists)
% -------------------------------------------------------------------------
out_filename = 'H3_04v3a_27_40';
% fprintf(1,'%s..',out_filename);
out_file0 = fullfile(cur_outpath,[out_filename,'.mat']);
if 0~=exist(out_file0,'file')
    return;
end

% -------------------------------------------------------------------------
% Panel initialization (peptide/PTM definitions, charges, m/z, RT seeds)
% -------------------------------------------------------------------------
His = init_histone(cur_outpath,out_filename);

% -------------------------------------------------------------------------
% Orchestration: relocation + extraction
% -------------------------------------------------------------------------
unitdiff = 1.0032;                  % C13-C12 mass difference (Da)
Mods = GetMods();                   % global PTM mass registry
[pep_rts,pep_intens,mono_isointens] = calculate_layout( ...
    MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,special);

% Persist quantitative results to disk
output_histone(cur_outpath,out_filename,His,pep_intens,pep_rts);

% Diagnostic plots (layout)
num_MS1 = size(MS1_index,1);
isorts = MS1_index(1:num_MS1,2);    % chronological RT vector (min)
draw_layout(cur_outpath,out_filename,His,pep_rts,pep_intens,isorts, ...
            mono_isointens,MS2_index,MS2_peaks,special);

% Optional PSM dump (DDA-like workflows)
if 1==special.nDAmode
    GetPSM(cur_outpath,out_filename,His,pep_rts,pep_intens,isorts, ...
           mono_isointens,MS1_index,MS1_peaks,MS2_index,ptol,unitdiff);
end

% ============================ NESTED FUNCTIONS ============================

function His = init_histone(cur_outpath,out_filename)
%%
% Build the panel for the v3 peptide sequence (KSAPSTGGVKKPHR).
% - mod_short: human-readable labels for each target PTM state
% - mod_type : compact encoding "idx,ptm;" where indices refer to residue
%              positions in the panel's internal mapping:
%              0 -> N-term (propionylation derivatization)
%              1 -> K27 (me1/2/3)
%              2 -> S28 (ph/ac/pr)
%             10 -> K36 (me1/2)
%             11 -> C-term derivatization
% - pep_ch   : practical charge states 2/3/4+
% - pep_mz   : theoretical m/z table computed by calculate_pepmz(His)
% - rt_ref   : seed RTs (min), later refined by relocate()/relocate2()
% - display  : mask to hide S28pr (channel 8) in figures by default

His.pep_seq = 'KSAPSTGGVKKPHR';

His.mod_short = {'S28ph';
    'K27me1S28ph';
    'K27me2S28ph';
    'K27me3S28ph';
    'K27me2S28phK36me2';
    'K27me3S28phK36me2';
    'S28ac';
    'S28pr';
    'K27me2S28ac';
    'K27me3S28ac';
    'K27me3S28acK36me1'};

His.mod_type = {'0,pr;1,pr;2,ph;10,pr;11,pr;';
    '0,pr;1,me1;2,ph;10,pr;11,pr;';
    '0,pr;1,me2;2,ph;10,pr;11,pr;';
    '0,pr;1,me3;2,ph;10,pr;11,pr;';
    '0,pr;1,me2;2,ph;10,me2;11,pr;';
    '0,pr;1,me3;2,ph;10,me2;11,pr;';
    '0,pr;1,pr;2,ac;10,pr;11,pr;';
    '0,pr;1,pr;2,pr;10,pr;11,pr;';
    '0,pr;1,me2;2,ac;10,pr;11,pr;';
    '0,pr;1,me3;2,ac;10,pr;11,pr;';
    '0,pr;1,me3;2,ac;10,me1;11,pr;'};

His.pep_ch = repmat([2 3 4],length(His.mod_type),1);

% Historical hard-coded tables are commented out above; prefer on-the-fly calc.
His.pep_mz = calculate_pepmz(His);

% Seed RTs (min) — copied from analogous panel; refined during relocation
His.rt_ref = [31.35
    33.20
    27.34
    27.31
    25
    25
    31.35
    33.20
    27.34
    27.31
    29.15];

% Hide S28pr by default in figures
His.display = ones(length(His.mod_type),1);
His.display(8) = 0;

His.outpath = cur_outpath;
His.outfile = out_filename;

% Normalize charge ordering so the "main" column (here 3+) is first
main_ch = His.pep_ch(1,2);
if main_ch~=His.pep_ch(1,1)
    [npep,ncharge] = size(His.pep_mz);
    new_ch = [main_ch,setdiff(His.pep_ch(1,:),main_ch)];
    x = zeros([1,ncharge]);
    for ino=1:ncharge
        x(ino) = find(His.pep_ch(1,:)==new_ch(ino));
    end
    tune = 1:npep;
    His.pep_mz(tune,:) = His.pep_mz(tune,x);
    His.pep_ch(tune,:) = His.pep_ch(tune,x);
end

function [pep_rts,pep_intens,mono_isointens] = calculate_layout( ...
    MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,Mods,His,special) %#ok
%%
% High-level controller:
%   1) Relocate RTs (debug/MS1/MS2 modes)
%   2) Extract RT apex + XIC intensities for all 11 target states via get_histone12
%
[npep,ncharge] = size(His.pep_mz);
num_MS1 = size(MS1_index,1);
pep_rts = zeros([npep,ncharge]);
pep_intens = zeros([npep,ncharge]);
mono_isointens = zeros([num_MS1,npep]);

% Relocation path
if 1==special.ndebug
    His = relocateD(MS1_index,MS1_peaks,ptol,unitdiff,His);
else
    if 2~=special.nDAmode
        His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His);
    else
        nhmass = special.nhmass;
        His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,nhmass);
    end
end

% Final extraction for all S28-bearing states
for hno=1:11
    [cur_rts,cur_intens,cur_mono_isointens] = get_histone12(MS1_index,MS1_peaks,ptol,unitdiff,His,hno);
    if cur_rts(1)>0
        pep_rts(hno,1:ncharge) = cur_rts;
        pep_intens(hno,1:ncharge) = cur_intens;
        mono_isointens(1:num_MS1,hno) = cur_mono_isointens;
    end
end

function His = relocate(MS1_index,MS1_peaks,ptol,unitdiff,His)
%%
% MS1-only relocation:
% - Read base RTs from H3_04v3_27_40.xls [rt] for KSAPSTGGVKKPHR.
% - Define windows for each S28-bearing target using those anchors.
% - Split me2 vs me3 (S28ph/ac) using find_pair() when both candidates exist.
%
delta = 0.3;
nsplit = 1;

ref_rts = get_KSAPSTGGVKKPHR_rt(His);
if 1==isempty(ref_rts)
    return;
end

% S28ph
hno = 1;
t1 = ref_rts(1)+delta;
t2 = ref_rts(1)+20;
[rts1,top1_rt1] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);
if 1==isempty(rts1), His.rt_ref(hno)=0; else, His.rt_ref(hno)=top1_rt1; end

% K27me1S28ph
hno = 2;
if His.rt_ref(1)>0,   t1 = His.rt_ref(1)+delta;
elseif ref_rts(3)>0,  t1 = ref_rts(3)+delta;
else,                 t1 = ref_rts(1)+delta; end
if His.rt_ref(1)>0,   t2 = His.rt_ref(1)+11;
elseif ref_rts(3)>0,  t2 = ref_rts(3)+18;
else,                 t2 = ref_rts(1)+25; end
[rts2,top1_rt2] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);
if 1==isempty(rts2), His.rt_ref(hno)=0; else, His.rt_ref(hno)=top1_rt2; end

% K27me2S28ph
hno = 3;
if 0==ref_rts(4), t1 = ref_rts(1)-35; else, t1 = ref_rts(4)+delta; end
if 0==ref_rts(9), t2 = ref_rts(1)-delta; else, t2 = ref_rts(9)-delta; end
[rts3,top1_rt3,inten_sum3,top1_inten_sum3] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

% K27me3S28ph
hno = 4;
if 0==ref_rts(4), t1 = ref_rts(1)-35; else, t1 = ref_rts(4)+delta; end
if 0==ref_rts(9), t2 = ref_rts(1)-delta; else, t2 = ref_rts(9)-delta; end
[rts4,top1_rt4,inten_sum4,top1_inten_sum4] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

[His.rt_ref(3),His.rt_ref(4)] = find_pair( ...
    rts3,top1_rt3,inten_sum3,top1_inten_sum3, ...
    rts4,top1_rt4,inten_sum4,top1_inten_sum4);

% K27me2S28phK36me2
hno = 5;
if 0==His.rt_ref(3), t1 = ref_rts(1)-45; else, t1 = His.rt_ref(3)-10; end
if 0==His.rt_ref(3), t2 = ref_rts(1)-10; else, t2 = His.rt_ref(3)-delta; end
[rts5,top1_rt5,inten_sum5,top1_inten_sum5] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

% K27me3S28phK36me2
hno = 6;
if 0==His.rt_ref(3), t1 = ref_rts(1)-45; else, t1 = His.rt_ref(3)-10; end
if 0==His.rt_ref(3), t2 = ref_rts(1)-10; else, t2 = His.rt_ref(3)-delta; end
[rts6,top1_rt6,inten_sum6,top1_inten_sum6] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

[His.rt_ref(5),His.rt_ref(6)] = find_pair( ...
    rts5,top1_rt5,inten_sum5,top1_inten_sum5, ...
    rts6,top1_rt6,inten_sum6,top1_inten_sum6);

% S28ac
hno = 7;
t1 = ref_rts(1)+delta;
if ref_rts(3)>0, t2 = ref_rts(3)+11; else, t2 = ref_rts(1)+18; end
[rts7,top1_rt7] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);
if 1==isempty(rts7), His.rt_ref(hno)=0; else, His.rt_ref(hno)=top1_rt7; end

% S28pr
hno = 8;
if His.rt_ref(7)>0,      t1 = His.rt_ref(7)+delta;
elseif ref_rts(3)>0,     t1 = ref_rts(3)+delta;
else,                    t1 = ref_rts(1)+delta; end
if His.rt_ref(7)>0,      t2 = His.rt_ref(7)+11;
elseif ref_rts(3)>0,     t2 = ref_rts(3)+18;
else,                    t2 = ref_rts(1)+25; end
[rts8,top1_rt8] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);
if 1==isempty(rts8), His.rt_ref(hno)=0; else, His.rt_ref(hno)=top1_rt8; end

% K27me2S28ac
hno = 9;
if ref_rts(4)>0, t1 = ref_rts(4)+delta; else, t1 = delta; end
t2 = ref_rts(1)-delta;
[rts9,top1_rt9,inten_sum9,top1_inten_sum9] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

% K27me3S28ac
hno = 10;
if ref_rts(4)>0, t1 = ref_rts(4)+delta; else, t1 = delta; end
t2 = ref_rts(1)-delta;
[rts10,top1_rt10,inten_sum10,top1_inten_sum10] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);

[His.rt_ref(9),His.rt_ref(10)] = find_pair( ...
    rts9,top1_rt9,inten_sum9,top1_inten_sum9, ...
    rts10,top1_rt10,inten_sum10,top1_inten_sum10);

% K27me3S28acK36me1
hno = 11;
if His.rt_ref(10)>0, t1 = His.rt_ref(10)+delta; else, t1 = delta; end
t2 = ref_rts(1)+delta;
[rts11,top1_rt11] = get_rts(MS1_index,MS1_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2);
if 1==isempty(rts11), His.rt_ref(hno)=0; else, His.rt_ref(hno)=top1_rt11; end

function His = relocate2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,nhmass)
%%
% MS2-assisted relocation:
% Same temporal windows as relocate(), but all searches use get_rts2() to
% enforce MS2 diagnostics (improved specificity in me2 vs me3 and K36 pairs).
%
delta = 0.3;
nsplit = 1;

ref_rts = get_KSAPSTGGVKKPHR_rt(His);
if 1==isempty(ref_rts)
    return;
end

% S28ph
hno = 1;
t1 = ref_rts(1)+delta;
t2 = ref_rts(1)+20;
[rts1,top1_rt1] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);
if 1==isempty(rts1), His.rt_ref(hno)=0; else, His.rt_ref(hno)=top1_rt1; end

% K27me1S28ph
hno = 2;
if His.rt_ref(1)>0,   t1 = His.rt_ref(1)+delta;
elseif ref_rts(3)>0,  t1 = ref_rts(3)+delta;
else,                 t1 = ref_rts(1)+delta; end
if His.rt_ref(1)>0,   t2 = His.rt_ref(1)+11;
elseif ref_rts(3)>0,  t2 = ref_rts(3)+18;
else,                 t2 = ref_rts(1)+25; end
[rts2,top1_rt2] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);
if 1==isempty(rts2), His.rt_ref(hno)=0; else, His.rt_ref(hno)=top1_rt2; end

% K27me2S28ph
hno = 3;
if 0==ref_rts(4), t1 = ref_rts(1)-35; else, t1 = ref_rts(4)+delta; end
if 0==ref_rts(9), t2 = ref_rts(1)-delta; else, t2 = ref_rts(9)-delta; end
[rts3,top1_rt3,inten_sum3,top1_inten_sum3] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

% K27me3S28ph
hno = 4;
if 0==ref_rts(4), t1 = ref_rts(1)-35; else, t1 = ref_rts(4)+delta; end
if 0==ref_rts(9), t2 = ref_rts(1)-delta; else, t2 = ref_rts(9)-delta; end
[rts4,top1_rt4,inten_sum4,top1_inten_sum4] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

[His.rt_ref(3),His.rt_ref(4)] = find_pair( ...
    rts3,top1_rt3,inten_sum3,top1_inten_sum3, ...
    rts4,top1_rt4,inten_sum4,top1_inten_sum4);

% K27me2S28phK36me2
hno = 5;
if 0==His.rt_ref(3), t1 = ref_rts(1)-45; else, t1 = His.rt_ref(3)-10; end
if 0==His.rt_ref(3), t2 = ref_rts(1)-10; else, t2 = His.rt_ref(3)-delta; end
[rts5,top1_rt5,inten_sum5,top1_inten_sum5] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

% K27me3S28phK36me2
hno = 6;
if 0==His.rt_ref(3), t1 = ref_rts(1)-45; else, t1 = His.rt_ref(3)-10; end
if 0==His.rt_ref(3), t2 = ref_rts(1)-10; else, t2 = His.rt_ref(3)-delta; end
[rts6,top1_rt6,inten_sum6,top1_inten_sum6] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

[His.rt_ref(5),His.rt_ref(6)] = find_pair( ...
    rts5,top1_rt5,inten_sum5,top1_inten_sum5, ...
    rts6,top1_rt6,inten_sum6,top1_inten_sum6);

% S28ac
hno = 7;
t1 = ref_rts(1)+delta;
if ref_rts(3)>0, t2 = ref_rts(3)+11; else, t2 = ref_rts(1)+18; end
[rts7,top1_rt7] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);
if 1==isempty(rts7), His.rt_ref(hno)=0; else, His.rt_ref(hno)=top1_rt7; end

% S28pr
hno = 8;
if His.rt_ref(7)>0,      t1 = His.rt_ref(7)+delta;
elseif ref_rts(3)>0,     t1 = ref_rts(3)+delta;
else,                    t1 = ref_rts(1)+delta; end
if His.rt_ref(7)>0,      t2 = His.rt_ref(7)+11;
elseif ref_rts(3)>0,     t2 = ref_rts(3)+18;
else,                    t2 = ref_rts(1)+25; end
[rts8,top1_rt8] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);
if 1==isempty(rts8), His.rt_ref(hno)=0; else, His.rt_ref(hno)=top1_rt8; end

% K27me2S28ac
hno = 9;
if ref_rts(4)>0, t1 = ref_rts(4)+delta; else, t1 = delta; end
t2 = ref_rts(1)-delta;
[rts9,top1_rt9,inten_sum9,top1_inten_sum9] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

% K27me3S28ac
hno = 10;
if ref_rts(4)>0, t1 = ref_rts(4)+delta; else, t1 = delta; end
t2 = ref_rts(1)-delta;
[rts10,top1_rt10,inten_sum10,top1_inten_sum10] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);

[His.rt_ref(9),His.rt_ref(10)] = find_pair( ...
    rts9,top1_rt9,inten_sum9,top1_inten_sum9, ...
    rts10,top1_rt10,inten_sum10,top1_inten_sum10);

% K27me3S28acK36me1
hno = 11;
if His.rt_ref(10)>0, t1 = His.rt_ref(10)+delta; else, t1 = delta; end
t2 = ref_rts(1)+delta;
[rts11,top1_rt11] = get_rts2(MS1_index,MS1_peaks,MS2_index,MS2_peaks,ptol,unitdiff,His,hno,nsplit,t1,t2,nhmass);
if 1==isempty(rts11), His.rt_ref(hno)=0; else, His.rt_ref(hno)=top1_rt11; end

function ref_rts = get_KSAPSTGGVKKPHR_rt(His)
%%
% Read base RTs for the v3 peptide KSAPSTGGVKKPHR from:
%   <His.outpath>/H3_04v3_27_40.xls
% The parser scans until the line "[rt]" is found, skips one header line,
% then parses subsequent TAB-separated lines and collects the 2nd column as RT.
out_file0 = fullfile(His.outpath,'H3_04v3_27_40.xls');
if 0~=exist(out_file0,'file')
    fp = fopen(out_file0,'r');
    str = fgetl(fp);
    while 0==feof(fp) && 0==strcmp(str,'[rt]')
        str = fgetl(fp);
    end
    str = fgetl(fp); % skip "peptide" header
    no = 0;
    ref_rts = [];
    while 0==feof(fp)
        str = fgetl(fp);
        p = strfind(str,'	');            % TABs
        c_rt = str2num( str(p(1)+1:p(2)-1) ); %#ok
        no = no + 1;
        ref_rts(no) = c_rt;               %#ok
    end
    fclose(fp);
else
    ref_rts = [];
end
