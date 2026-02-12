function draw_layout(cur_outpath,out_filename,His,pep_rts,pep_intens,isorts,mono_isointens,MS2_index,MS2_peaks,special)
%%

npep = size(His.pep_mz,1);
ix = find(pep_rts(1:npep,1)>4 & His.display==1);
if 1==isempty(ix)
    return;
end;

% localmax_rt, localmax_inten, terminus
nplot = length(ix);
localmax_rt = zeros([nplot,1]);
localmax_inten = zeros([nplot,1]);
terminus = zeros([nplot,2]);
for ino=1:nplot
    cno = ix(ino);
    p = find(isorts<=pep_rts(cno,1));
    c_ms1pos = p(end);
    c_mono_isointens = mono_isointens(:,cno);
    if pep_intens(cno,1)>0
        [nt,nb] = GetTopBottom11(c_mono_isointens);%#ok
        [localmax_rt(ino),localmax_inten(ino),IX] = GetLocal(c_ms1pos,isorts,c_mono_isointens,nb);
        if 1==isempty(IX)
            terminus(ino,1:2) = [c_ms1pos c_ms1pos];
        else
            terminus(ino,1:2) = [IX(1) IX(end)];
        end;
        localmax_rt(ino) = pep_rts(cno,1);
    else
        terminus(ino,1:2) = [c_ms1pos c_ms1pos];
        localmax_rt(ino) = pep_rts(cno,1);
    end;
end;
maxinten = max(localmax_inten);%#ok
st = max([isorts(min(terminus(1:nplot,1)))-5 1]);
tm = isorts(max(terminus(1:nplot,2)))+25;

out_file1 = fullfile(fileparts(cur_outpath),[out_filename,'.pdf']);
warning off all;
set(gcf,'visible','off');
% set(gcf,'position',[0 0 1000 1000]);
if 1==strcmp(out_filename(1:2),'HH')
    out_filename = out_filename(2:end);
end;
p = strfind(out_filename,'_');
cur_title = [out_filename(1:p(1)-1),' ',out_filename(p(2)+1:p(3)-1),'-',out_filename(p(3)+1:end),' ',His.pep_seq,' +',num2str(His.pep_ch(1,1)),' ions'];
nhmass = special.nhmass;

Mods = GetMods();
colors = {'k','r','g','b','c','m'};

for ino=1:nplot
    cno = ix(ino);
    subplot(nplot,1,ino);
    % XIC
    plot(isorts,mono_isointens(:,cno),'color','b','linewidth',1);
    set(gca,'xtick',[],'ytick',[]);
    hold on;
    xlim([st tm]);
    %ylim([0 1.05*maxinten]);

    % --- Robustez (tercer error): mapear [st, tm] a índices de isorts sin fallar ---
    tmp_maxinten = maxinten; % valor por defecto seguro si no podemos recalcular

    p1 = find(isorts <= st);
    p2 = find(isorts <= tm);

    if isempty(p1) || isempty(p2)
        % Si st < min(isorts) o tm NaN/fuera de rango, evitamos p1(end)/p2(end)
        fprintf('draw_layout: rango vacío st=%g tm=%g RT=[%g,%g]\n', st, tm, min(isorts), max(isorts));
    else
        i1 = p1(end);
        i2 = p2(end);
        if i2 >= i1
            IX = i1:i2;
            tmp_maxinten = max(mono_isointens(IX,cno));
            if tmp_maxinten > 0
                ylim([0 1.05*tmp_maxinten]);
            end
        else
            fprintf('draw_layout: rango invertido st=%g tm=%g i1=%d i2=%d\n', st, tm, i1, i2);
        end
    end
    % --- fin robustez ---

    % localmax
    plot(localmax_rt(ino),localmax_inten(ino),'color','m','linestyle','-','linewidth',1);
    cur_txt = [His.mod_short{cno},'(',num2str(His.pep_mz(cno,1),'%.4f'),',+',num2str(His.pep_ch(cno,1)),'), ',num2str(localmax_rt(ino),'%.2f'),', ',num2str(localmax_inten(ino),'%.2e')];
    %text(localmax_rt(ino),localmax_inten(ino)+0.05*maxinten,cur_txt,'color','r','fontsize',7);
    if pep_intens(cno,1)>0
        text(localmax_rt(ino),1.05*tmp_maxinten,cur_txt,'color','r','fontsize',7);
    else
        text(localmax_rt(ino),localmax_inten(ino),cur_txt,'color','r','fontsize',7);
    end;

    % boundary
    plot([isorts(terminus(ino,1)) isorts(terminus(ino,1))],[0 1],'color','m','linestyle','-','linewidth',1);
    plot([isorts(terminus(ino,2)) isorts(terminus(ino,2))],[0 1],'color','m','linestyle','-','linewidth',1);
    if 1==ino
        title(cur_title);
    end;

    % DIA
    if 2==special.nDAmode
        % match MS2
        rt1 = isorts(terminus(ino,1))-0.001;
        rt2 = isorts(terminus(ino,2))+0.001;
        [ms2pos,ms2rts,ms2intens,posn,posc,ActiveType] = MatchMS2(MS2_index,MS2_peaks,Mods,His,cno,rt1,rt2,nhmass);
        if 1==isempty(ms2pos)
            continue;
        end;

        % fragment ions
        if 1==strcmp(ActiveType,'CID')
            strn = 'b';
            strc = 'y';
        else
            strn = 'c';
            strc = 'z';
        end;
        new_maxinten = max(max(ms2intens));
        if new_maxinten>0
            fold = (tmp_maxinten/new_maxinten)/(1+length(posn));
        else
            fold = 1/(1+length(posn));
        end;
        for kno=1:length(posn)
            plot(ms2rts(ms2pos)-4,fold*ms2intens(ms2pos,kno)+(kno-1)*fold*new_maxinten,'color',colors{mod(kno,6)+1},'linestyle','-','linewidth',0.5);
            text(ms2rts(ms2pos(end))-4,fold*ms2intens(ms2pos(end),kno)+(kno-1)*fold*new_maxinten,[strn,num2str(posn(kno))],'color',colors{mod(kno,6)+1},'fontsize',7);
        end;
        for kno=1:length(posc)
            qno = kno+length(posn);
            plot(ms2rts(ms2pos)+3,fold*ms2intens(ms2pos,qno)+(kno-1)*fold*new_maxinten,'color',colors{mod(kno,6)+1},'linestyle','-','linewidth',0.5);
            text(ms2rts(ms2pos(end))+3,fold*ms2intens(ms2pos(end),qno)+(kno-1)*fold*new_maxinten,[strc,num2str(posc(kno))],'color',colors{mod(kno,6)+1},'fontsize',7);
        end;
    end;
end;
set(gca,'xtickMode', 'auto');
xlabel('Time (min)');
ylabel('Abundance');
print('-dpdf',out_file1);
close();

function [ms2pos,ms2rts,ms2intens,posn,posc,ActiveType] = MatchMS2(MS2_index,MS2_peaks,Mods,His,hno,rt1,rt2,nhmass)
%%

% get precursors in MS1 profile
num_MS2 = size(MS2_index,1);
c_mz = His.pep_mz(hno,1);
premzs = unique(MS2_index(:,4));
[tmp,ii] = min( abs(premzs-c_mz) );%#ok
target = premzs(ii);
flag = zeros([num_MS2,1]);
p = find( MS2_index(:,2)>=rt1 );
pp = find( MS2_index(:,2)<=rt2 );
if 1==isempty(p) || 1==isempty(pp) || p(1)>pp(end)
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
for i=i1:i2
    cen_mz = MS2_index(i,4);
    if 0==cen_mz-target
        flag(i) = 1;
    end;
end;
ms2pos = find(flag==1);
if 1==isempty(ms2pos)
    ms2pos = i1:i2;
end;

% check MS2
ms2rts = MS2_index(:,2);

instruments = MS2_index(ms2pos,6);% MS2dirs = {'CIDIT','CIDFT','ETDIT','ETDFT','HCDIT','HCDFT'};
% if 1==length(unique(instruments))
    % ActiveType, tol
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

    % K1,K2
    if 1==nhmass
        [K1,posn,posc] = get_key_ions1H(His,hno,Mods,ActiveType);
    else
        [K1,posn,posc] = get_key_ions1(His,hno,Mods,ActiveType);
    end;
% end;

index = [1;MS2_index(1:num_MS2,7)];
ms2intens = zeros([num_MS2,length(K1)]);
for i=1:length(ms2pos)
    cno = ms2pos(i);
    for pno = cno%cno-1:cno+1
        if pno<1 || pno>num_MS2
            continue;
        end;
        if 1<length(unique(instruments))
            % ActiveType, tol
            c_instrument = MS2_index(pno,6);% MS2dirs = {'CIDIT','CIDFT','ETDIT','ETDFT','HCDIT','HCDFT'};
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

        if 1<length(unique(instruments))
            % K1,K2
            if 1==nhmass
                [K1,posn,posc] = get_key_ions1H(His,hno,Mods,ActiveType);
            else
                [K1,posn,posc] = get_key_ions1(His,hno,Mods,ActiveType);
            end;
        end;

        % mz, inten
        IX = index(pno):index(pno+1)-1;
        mz = MS2_peaks(IX,1);
        inten = MS2_peaks(IX,2);

        % match key ions
        for j=1:length(K1)
            ix1 = find(abs(mz-K1(j))<=tol);
            [tmp,x1] = min(abs(mz(ix1)-K1(j)));%#ok
            %[tmp,x1] = max(inten(ix1));
            if 0==isempty(ix1) && ms2intens(cno,j)<inten(ix1(x1))
                ms2intens(cno,j) = inten(ix1(x1));
            end;
        end;
    end;
end;
