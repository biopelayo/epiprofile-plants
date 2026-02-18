function success = GetMS1ScanNo(ms1_fullfile)
%% GetMS1ScanNo

success = 0;

% check the MAT file
[datapath,dataname] = fileparts(ms1_fullfile);
MS1_scanfile = fullfile(datapath,[dataname,'_MS1scans.mat']);
MS1_peakfile = fullfile(datapath,[dataname,'_MS1peaks.mat']);
if 0~=exist(MS1_scanfile,'file') && 0~=exist(MS1_peakfile,'file')
    %{%for old data do not have 'MS1Type'
    load(MS1_scanfile);
    if 0==exist('MS1Type')%#ok
        MS1Type = 'FTMS';%#ok
        save(MS1_scanfile,'MS1_index','MS1Type');
    end;
    %}
    success = 1;
    return;
end;

% check the TXT file
if 0==exist(ms1_fullfile,'file')
    disp([ms1_fullfile,': does not exist!']);
    return;
end;

% open the TXT file
fid = fopen(ms1_fullfile,'r');
if -1==fid
    disp([ms1_fullfile,': can not open!']);
    return;
end;

%% init
% initialize the MS1 info
%ptol = 20;
maxpeaknum = 1e4;% the max peak number on a MS1 scan
maxMS1num = 1.5e5;% initial MS1 scan number
totalpeaknum = 4e7;% the init total peak number on MS1 scans
MS1_index = zeros([maxMS1num,4]);% MS1 scan, MS1 rt, MS1 peak num, baseline
MS1_peaks = zeros([totalpeaknum,2]);% m/z and intensity on MS1
fno = 0;% real MS1 scan number
pkno = 0;% real total peak number

% get the keywords
keyword0 = 'H	DataType';
keyword1 = 'S';% the keyword to start record
keyword2 = 'I	RetTime';% the keyword to start record
keyword3 = 'I	InstrumentType';
len0 = length(keyword0);
len1 = length(keyword1);
len2 = length(keyword2);
len3 = length(keyword3);

%% get the MS1 info
% get the datatype
str=fgets(fid);
while feof(fid)==0 && 0==strcmp( str(1:len0),keyword0 )
    str=fgets(fid);
end;
MS1_datamode = str(len0+2);

if 1==strcmp('P',MS1_datamode)
    fprintf(1,'MS1 is profile mode, convert to centroid mode first!\n');
    fclose(fid);
    return;
end;

% for progress
ct_prt = 0;
fprintf(1,'MS1 scans: ');

% start to process
str = fgets(fid);
while 0==feof(fid)
    if 1==strcmp( str(1:len1),keyword1 )
        % progress
        fno = fno + 1;
        fprintf(repmat('\b',[1,ct_prt]));
        ct_prt = fprintf('%i',fno);

        % 1.get the MS1 info
        % MS1 scan
        scan_no = str2num(str(len1+2:end));%#ok
        scan_no = scan_no(1);

        % RT
        str = fgets(fid);
        while feof(fid)==0 && 0==strcmp( str(1:len2),keyword2 )
            str=fgets(fid);
        end;
        rt_no = str2double(str(len2+2:end));

        % IonInjectionTime
        str = fgets(fid);%#ok

        % InstrumentType
        str = fgets(fid);
		%
        if 1==fno
            MS1Type = str(len3+2:len3+5);%#ok
%             if 1==strcmp(MS1Type,'ITMS')
%                 fclose(fid);
%                 disp([ms1_fullfile,': low resolution MS1!']);
%                 return;
%             end;
        end;
		%}

        % 2.read the MS1 data
        mz = zeros([1,maxpeaknum]);
        inten = zeros([1,maxpeaknum]);
        pnum = 0;
        while feof(fid)==0 && 0==strcmp( str(1:len1),keyword1 )
            if ~('0'<=str(1) && str(1)<='9')
                str = fgets(fid);
                continue;
            end;
            pnum = pnum + 1;
            tmp = sscanf(str,'%f',2);
            mz(pnum) = tmp(1);
            inten(pnum) = tmp(2);
            str = fgets(fid);
        end;
        if 1==feof(fid)
            pnum = pnum + 1;
            tmp = sscanf(str,'%f',2);
            mz(pnum) = tmp(1);
            inten(pnum) = tmp(2);
        end;
        IX = find(inten>0);% first read, then filter
        mz = mz(IX);
        inten = inten(IX);

        % 3.judge whether to centroid MS1 scan

        % 4.save the MS1 info and peaks
        if 1==isempty(mz) || mz(end)-mz(1)<10
            fno = fno - 1;
            continue;
        end;
        if length(inten)>100
            baseline = GetBaseline(inten);
            IX = find(inten>=baseline);
            mz = mz(IX);
            inten = inten(IX);
        else
            baseline = 0;
        end;
        npk = length(mz);

        MS1_index(fno,1:4) = [scan_no rt_no npk baseline];
        MS1_peaks(pkno+1:pkno+npk,1) = mz;
        MS1_peaks(pkno+1:pkno+npk,2) = inten;
        pkno = pkno + npk;
    else
        str = fgets(fid);
    end;
end;
fclose(fid);% close the TXT file
fprintf(repmat('\b',[1,ct_prt]));
fprintf('%i',fno);
fprintf(1,'\n');

%% save the MS1 info
% filter the empty values
if fno<maxMS1num
    IX = 1:fno;
    MS1_index = MS1_index(IX,:);
end;
tmp = MS1_index(1:fno,3);
MS1_index(1:fno,3) = cumsum(tmp) + 1;
if MS1_index(fno,2)>1000% 20*60=1200
    MS1_index(1:fno,2) = MS1_index(1:fno,2)/60;%#ok
end;

if pkno<totalpeaknum
    IX = 1:pkno;
    MS1_peaks = MS1_peaks(IX,:);%#ok
end;

% save the results
save(MS1_scanfile,'MS1_index','MS1Type');
save(MS1_peakfile,'MS1_peaks');

success = 1;

function baseline = GetBaseline(inten)
%% GetBaseline

loginten = log10(inten);
t = min(loginten):0.08:max(loginten);
[n,xout] = hist(loginten,t);

[tmp,idx] = max(n);%#ok
baseline = 10^xout(idx);