function draw_layout(cur_outpath,out_filename,His,pep_rts,pep_intens,isorts,mono_isointens,MS2_index,MS2_peaks,special)
%% ========================================================================
% draw_layout
% -------------------------------------------------------------------------
% PURPOSE (plain English, non-proteomics audience):
%   Produce a per-peptide visualization (as a multi-panel PDF) of MS1
%   chromatographic traces (XICs) for all modified forms of a peptide, and
%   optionally overlay diagnostic MS2 fragment evidence (DIA-like mode).
%
% WHAT IT DOES:
%   1) Chooses which peptide rows to plot (filters by RT and a display flag).
%   2) For each selected row:
%        - Plots the MS1 extracted ion chromatogram (mono-isotopic intensity).
%        - Marks a local maximum and the integration boundaries in time.
%        - Writes an informative text label (modification, m/z, charge, RT, intensity).
%        - Optionally overlays MS2 fragment intensities aligned to time.
%   3) Saves a single PDF with N subplots (one per selected row).
%
% INPUTS (high-level meanings):
%   cur_outpath     : folder path where detail files live (used to place the PDF).
%   out_filename    : logical peptide/window name used for PDF title & output name.
%   His             : struct containing peptide metadata for this window.
%                     Required fields (as used here):
%                       - pep_mz   : [nRows x 1] precursor m/z per row
%                       - pep_ch   : [nRows x 1] charge state per row
%                       - mod_short: {nRows x 1} short label of the mod per row
%                       - pep_seq  : (char) peptide sequence for the window
%                       - display  : [nRows x 1] flag; rows with 1 are considered plottable
%   pep_rts         : [nRows x 1] nominal RT per row (minutes).
%   pep_intens      : [nRows x 1] a scalar intensity per row (used to gate labels).
%   isorts          : [nTimePoints x 1] time axis (minutes) for MS1 profile.
%   mono_isointens  : [nTimePoints x nRows] mono-isotopic MS1 intensities per row.
%   MS2_index       : index table for MS2 scans (instrument metadata).
%                     Columns used here:
%                       (2) RT (minutes)
%                       (4) precursor m/z
%                       (6) instrument code (1..6)
%                       (7) cumulative/row index for peaks in MS2_peaks (see below)
%   MS2_peaks       : concatenated MS2 m/z-intensity pairs for all scans.
%                     For a given scan 's', its rows are slice
%                     MS2_peaks(index(s):index(s+1)-1, :).
%   special         : struct with mode/chemistry flags:
%                       - nDAmode : if ==2, overlay MS2 fragments (DIA-like)
%                       - nhmass  : use H-mass variant when choosing key ions
%
% OUTPUT:
%   Writes a PDF "<parent_of_cur_outpath>/<out_filename>.pdf" with the plot.
%
% EXTERNAL HELPERS (not included here):
%   GetTopBottom11, GetLocal, GetMods, get_key_ions1, get_key_ions1H
%   and the local nested function MatchMS2 (defined below).
%
% IMPORTANT VISUAL NOTES:
%   - X axis is time (minutes), Y axis is intensity (arbitrary units).
%   - Every subplot corresponds to one modified form (row) of the peptide.
%   - Magenta vertical lines delineate the integration bounds for each row.
%   - When DIA overlay is active, colored small traces/labels show fragment
%     ions (n-/c-series depending on instrument: CID→b/y, ETD→c/z).
% ========================================================================

% ----------------------------
% Basic selection of peptide rows to plot:
%   npep  = number of rows (modified forms) in this peptide window.
%   ix    = indices (rows) where: RT > 4 minutes AND His.display == 1.
%           (RT gate removes early noise; display gate respects prior QC.)
npep = size(His.pep_mz,1);
ix = find(pep_rts(1:npep,1)>4 & His.display==1);
if 1==isempty(ix)
    % Nothing to plot for this window; exit silently.
    return;
end;

% ----------------------------
% Preallocate containers for per-row local maxima and boundaries:
%   localmax_rt/inten : a representative apex in time & intensity.
%   terminus          : [startIndex endIndex] in 'isorts' defining the
%                       integration window for this row.
nplot = length(ix);
localmax_rt = zeros([nplot,1]);
localmax_inten = zeros([nplot,1]);
terminus = zeros([nplot,2]);

% For each selected row:
for ino=1:nplot
    cno = ix(ino);  % current row index in this window

    % Find the closest time index in 'isorts' that is <= the row RT.
    % c_ms1pos will seed local searching around the nominal RT.
    p = find(isorts<=pep_rts(cno,1));
    c_ms1pos = p(end);

    % Column slice of the mono-isotopic XIC for this row.
    c_mono_isointens = mono_isointens(:,cno);

    % If this row has positive scalar intensity, try to find local apex and bounds.
    if pep_intens(cno,1)>0
        [nt,nb] = GetTopBottom11(c_mono_isointens); %#ok  % returns top/bottom guides (not reused)
        [localmax_rt(ino),localmax_inten(ino),IX] = GetLocal(c_ms1pos,isorts,c_mono_isointens,nb);

        % If GetLocal returned no boundary indices, fall back to a degenerate bound.
        if 1==isempty(IX)
            terminus(ino,1:2) = [c_ms1pos c_ms1pos];
        else
            terminus(ino,1:2) = [IX(1) IX(end)];
        end;

        % IMPORTANT: Regardless of the apex returned by GetLocal,
        % lock the displayed apex time to the nominal peptide RT.
        localmax_rt(ino) = pep_rts(cno,1);
    else
        % If the row has no scalar intensity, mark a degenerate window at the seed.
        terminus(ino,1:2) = [c_ms1pos c_ms1pos];
        localmax_rt(ino) = pep_rts(cno,1);
    end;
end;

% Global max (kept for potential scaling; not directly used below).
maxinten = max(localmax_inten); %#ok

% Global X limits for the figure:
%   st = left bound (a bit before the earliest left terminus)
%   tm = right bound (a bit after the latest right terminus)
st = max([isorts(min(terminus(1:nplot,1)))-5 1]);
tm = isorts(max(terminus(1:nplot,2)))+25;

% ----------------------------
% Output PDF path: one level above 'cur_outpath', named after 'out_filename'.
out_file1 = fullfile(fileparts(cur_outpath),[out_filename,'.pdf']);

% Figure setup (headless):
warning off all;
set(gcf,'visible','off');
% set(gcf,'position',[0 0 1000 1000]);  % (optional: figure size)

% Normalize an 'HH' prefix in out_filename (historical naming quirk).
if 1==strcmp(out_filename(1:2),'HH')
    out_filename = out_filename(2:end);
end;

% Build a readable title from filename tokens + peptide sequence + charge.
p = strfind(out_filename,'_');
cur_title = [out_filename(1:p(1)-1),' ',out_filename(p(2)+1:p(3)-1),'-',out_filename(p(3)+1:end),' ',His.pep_seq,' +',num2str(His.pep_ch(1,1)),' ions'];

% Chemistry/switches carried in 'special'.
nhmass = special.nhmass;

% Load modification catalog for MS2 matching (external helper).
Mods = GetMods();

% Color palette for fragment overlays (cycles over 6 colors).
colors = {'k','r','g','b','c','m'};

% ----------------------------
% MAIN PLOTTING LOOP: one subplot per selected row (modified form).
for ino=1:nplot
    cno = ix(ino);
    subplot(nplot,1,ino);

    % === 1) MS1 XIC trace ===
    plot(isorts,mono_isointens(:,cno),'color','b','linewidth',1);
    set(gca,'xtick',[],'ytick',[]);
    hold on;
    xlim([st tm]);

    % Adapt Y limits to the local window [st, tm] for better readability.
    p1 = find(isorts<=st);
    p2 = find(isorts<=tm);
    IX = p1(end):p2(end);
    tmp_maxinten = max(mono_isointens(IX,cno));
    if tmp_maxinten>0
        ylim([0 1.05*tmp_maxinten]);
    end;

    % === 2) Mark the "local maximum" and annotate with a text label ===
    plot(localmax_rt(ino),localmax_inten(ino),'color','m','linestyle','-','linewidth',1);
    cur_txt = [His.mod_short{cno},'(',num2str(His.pep_mz(cno,1),'%.4f'),',+',num2str(His.pep_ch(cno,1)),'), ',num2str(localmax_rt(ino),'%.2f'),', ',num2str(localmax_inten(ino),'%.2e')];

    % Place the text above the trace if there is signal; otherwise at apex point.
    if pep_intens(cno,1)>0
        text(localmax_rt(ino),1.05*tmp_maxinten,cur_txt,'color','r','fontsize',7);
    else
        text(localmax_rt(ino),localmax_inten(ino),cur_txt,'color','r','fontsize',7);
    end;

    % === 3) Draw vertical magenta lines at integration boundaries ===
    plot([isorts(terminus(ino,1)) isorts(terminus(ino,1))],[0 1],'color','m','linestyle','-','linewidth',1);
    plot([isorts(terminus(ino,2)) isorts(terminus(ino,2))],[0 1],'color','m','linestyle','-','linewidth',1);

    % Title only on the first subplot.
    if 1==ino
        title(cur_title);
    end;

    % === 4) Optional DIA-like MS2 overlay (if special.nDAmode == 2) ===
    if 2==special.nDAmode
        % Restrict MS2 matching to the same integration window with slight padding.
        rt1 = isorts(terminus(ino,1))-0.001;
        rt2 = isorts(terminus(ino,2))+0.001;

        % Try to match fragment ions in MS2 within [rt1, rt2].
        % Returns:
        %   ms2pos    : scan indices to draw
        %   ms2rts    : RT vector for all scans
        %   ms2intens : [nScans x nFragments] intensity table
        %   posn/posc : fragment indices for N- and C-terminal series
        %   ActiveType: 'CID' → b/y series; 'ETD' → c/z series
        [ms2pos,ms2rts,ms2intens,posn,posc,ActiveType] = MatchMS2(MS2_index,MS2_peaks,Mods,His,cno,rt1,rt2,nhmass);
        if 1==isempty(ms2pos)
            % No MS2 overlay to draw in this window.
            continue;
        end;

        % Choose fragment letters according to activation type.
        if 1==strcmp(ActiveType,'CID')
            strn = 'b';  % N-terminal series
            strc = 'y';  % C-terminal series
        else
            strn = 'c';
            strc = 'z';
        end;

        % Compute a vertical scaling 'fold' so fragment tracks fit the panel.
        new_maxinten = max(max(ms2intens));
        if new_maxinten>0
            fold = (tmp_maxinten/new_maxinten)/(1+length(posn));
        else
            fold = 1/(1+length(posn));
        end;

        % Draw N-series fragments shifted slightly left in time (−4 min).
        for kno=1:length(posn)
            plot(ms2rts(ms2pos)-4, ...
                 fold*ms2intens(ms2pos,kno)+(kno-1)*fold*new_maxinten, ...
                 'color',colors{mod(kno,6)+1},'linestyle','-','linewidth',0.5);
            text(ms2rts(ms2pos(end))-4, ...
                 fold*ms2intens(ms2pos(end),kno)+(kno-1)*fold*new_maxinten, ...
                 [strn,num2str(posn(kno))], ...
                 'color',colors{mod(kno,6)+1},'fontsize',7);
        end;

        % Draw C-series fragments shifted slightly right in time (+3 min).
        for kno=1:length(posc)
            qno = kno+length(posn);
            plot(ms2rts(ms2pos)+3, ...
                 fold*ms2intens(ms2pos,qno)+(kno-1)*fold*new_maxinten, ...
                 'color',colors{mod(kno,6)+1},'linestyle','-','linewidth',0.5);
            text(ms2rts(ms2pos(end))+3, ...
                 fold*ms2intens(ms2pos(end),qno)+(kno-1)*fold*new_maxinten, ...
                 [strc,num2str(posc(kno))], ...
                 'color',colors{mod(kno,6)+1},'fontsize',7);
        end;
    end;
end;

% Restore a visible X-axis on the last axes and add global labels.
set(gca,'xtickMode', 'auto');
xlabel('Time (min)');
ylabel('Abundance');

% Export to PDF and close the figure.
print('-dpdf',out_file1);
close();


% ========================================================================
% NESTED HELPER: MatchMS2
% ------------------------------------------------------------------------
% PURPOSE:
%   Given a time window and a target precursor m/z, retrieve MS2 scans,
%   decide activation type, select key fragment ions, and extract their
%   matched intensities per scan to be overlaid in the plot.
% ========================================================================
function [ms2pos,ms2rts,ms2intens,posn,posc,ActiveType] = MatchMS2(MS2_index,MS2_peaks,Mods,His,hno,rt1,rt2,nhmass)
%%

% ----------------------------
% Determine which MS2 scans in the window correspond to the target precursor:
num_MS2 = size(MS2_index,1);
c_mz = His.pep_mz(hno,1);

% Unique precursor m/z values present in MS2_index:
premzs = unique(MS2_index(:,4));

% Choose the closest present precursor to the current target m/z:
[tmp,ii] = min( abs(premzs-c_mz) ); %#ok
target = premzs(ii);

% Filter scans by time within [rt1, rt2].
flag = zeros([num_MS2,1]);
p = find( MS2_index(:,2)>=rt1 );
pp = find( MS2_index(:,2)<=rt2 );
if 1==isempty(p) || 1==isempty(pp) || p(1)>pp(end)
    % No MS2 scans in that window → early exit with empties.
    ms2pos = [];
    ms2rts = [];
    ms2intens = [];
    posn = [];
    posc = [];
    ActiveType = [];
    return;
end;
i1 = p(1);
i2 = pp(end);

% Mark scans whose reported precursor m/z equals the chosen 'target'.
for i=i1:i2
    cen_mz = MS2_index(i,4);
    if 0==cen_mz-target
        flag(i) = 1;
    end;
end;

% Positions of "matching" scans:
ms2pos = find(flag==1);
if 1==isempty(ms2pos)
    % If no exact match, fallback to all scans in the window (lenient).
    ms2pos = i1:i2;
end;

% ----------------------------
% Prepare outputs and instrument-type logic:
ms2rts = MS2_index(:,2);

% instruments: numeric codes (1..6). Mapping (comment from original code):
%   {'CIDIT','CIDFT','ETDIT','ETDFT','HCDIT','HCDFT'}
instruments = MS2_index(ms2pos,6);

% Decide activation type and tolerance from the FIRST instrument code:
%   CID/HCD → use 'CID' branch (b/y series)
%   ETD     → use 'ETD' branch (c/z series)
%   Odd codes → "ion trap" style → wider tol (0.4)
%   Even codes → "FT" style       → tighter tol (0.02)
% NOTE: If scans are heterogeneous, we re-decide inside the loop below.
c_instrument = instruments(1);
if 3==c_instrument || 4==c_instrument
    ActiveType = 'ETD';
else
    ActiveType = 'CID';
end;
if 1==mod(c_instrument,2)
    tol = 0.4;
else
    tol = 0.02;
end;

% Choose key fragment m/z list (K1) and their positions for labeling.
% nhmass==1 → use the H-mass variant helper.
if 1==nhmass
    [K1,posn,posc] = get_key_ions1H(His,hno,Mods,ActiveType);
else
    [K1,posn,posc] = get_key_ions1(His,hno,Mods,ActiveType);
end;

% ----------------------------
% Build an index (row starts per scan) to slice MS2_peaks per scan:
% 'index' is 1-based, with an extra sentinel at the end via MS2_index(:,7).
index = [1;MS2_index(1:num_MS2,7)];

% Preallocate matrix of matched fragment intensities:
%   rows   → scan index (global)
%   cols   → fragment m/z entries from K1 (N- then C-series)
ms2intens = zeros([num_MS2,length(K1)]);

% For each selected scan position:
for i=1:length(ms2pos)
    cno = ms2pos(i);

    % Consider the "central" scan only (as per original code).
    for pno = cno % (legacy: could be cno-1:cno+1)
        if pno<1 || pno>num_MS2
            continue;
        end;

        % If instruments are heterogeneous across scans, re-compute type/tol:
        if 1<length(unique(instruments))
            c_instrument = MS2_index(pno,6);
            if 3==c_instrument || 4==c_instrument
                ActiveType = 'ETD';
            else
                ActiveType = 'CID';
            end;
            if 1==mod(c_instrument,2)
                tol = 0.4;
            else
                tol = 0.02;
            end;
        end;

        % If heterogeneous instruments, re-select K1 accordingly as well.
        if 1<length(unique(instruments))
            if 1==nhmass
                [K1,posn,posc] = get_key_ions1H(His,hno,Mods,ActiveType);
            else
                [K1,posn,posc] = get_key_ions1(His,hno,Mods,ActiveType);
            end;
        end;

        % Retrieve the m/z-intensity pairs for this MS2 scan 'pno'.
        IX = index(pno):index(pno+1)-1;
        mz = MS2_peaks(IX,1);
        inten = MS2_peaks(IX,2);

        % Match each target fragment m/z in K1 within absolute tolerance 'tol'.
        for j=1:length(K1)
            ix1 = find(abs(mz-K1(j))<=tol);
            [tmp,x1] = min(abs(mz(ix1)-K1(j))); %#ok
            % Optionally, could choose max intensity instead of nearest m/z:
            % [~,x1] = max(inten(ix1));
            if 0==isempty(ix1) && ms2intens(cno,j)<inten(ix1(x1))
                ms2intens(cno,j) = inten(ix1(x1));
            end;
        end;
    end;
end;
