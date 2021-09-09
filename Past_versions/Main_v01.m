%% batch_res_v14 - batch analysis of radial pulse data using kreservoir_v14
%% Copyright 2020 Alun Hughes & Kim Parker
%% Modified by Lin for Project evluating effectiveness of pressure-only WIA against normal WIA.
% This software is distributed under under the terms of the GNU General Public License
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% http://www.gnu.org/licenses/gpl.html

%% Versions
%  v0.1 (First version working for virtual population with pulse and signal input)


%%%%%%%%%%%%%%%%
%% m files required to be in directory
% fitres_v6.m
% kreservoir_v14.m
% hrv_v1.m
%%%%%%%%%%%%%%%%
%% Constants
clear all; clc;

sampling_rate = 1000;  % For our sample interval of 1ms.

kres_v='v14';          % Version tracking for reservoir fitting
headernumber=47;       % headers for columns of results (see end)
mmHgPa = 133;          % P conversion for WIA
uconst=1;              % empirical constant to convert normalized
% velocity to m/s based on Hughes et al.
% Front. Physiol. DOI: 10.3389/fphys.2020.00550
Npoly=3;               % Order of polynomial fit for sgolay
Frame=9;               % Window length for sgolay
version='0.1';         % Version of Main
%%%%%%%%%%%%%%%%
%% Select files
% default folder as per manual
folder_name ='C:\Spdata\';
% check that C:\Spdata exists and if not allows new folder to be chosen
if ~exist('C:\Spdata', 'dir')
    answer = questdlg('C:\Spdata doesnt exist. Would you like to choose another folder?', ...
        'Sphygmocor Data Folder','Yes', 'No [end]','Yes');
    % Handle response
    switch answer
        case 'Yes'
            folder_name = uigetdir;
        case 'No [end]'             % end if no folder identified
            return
    end
end
% read files
file_lists=dir(fullfile(folder_name, '*.txt'));
no_of_files=length(file_lists);
% add an error trap here if no files in folder
if no_of_files==0
    f = errordlg('No data files to analyse in folder','File error');
    return
end
% set record number to 1 and extract filename
record_no=1;
% preallocate cell array
proc_var=cell(no_of_files,headernumber);

%% Reservoir analysis for each filename
for file_number=1:no_of_files
    % refresh filename
    filename=file_lists(record_no).name;
    % read BP data from Sphygmocor file using textread as it's easier
    %       [periph_signal,central_signal,periph_pulse,central_pulse,...
    %           flow_waveform,forward_pulse,reflected_pulse] = textread...
    %           ([folder_name filename],'%f%f%f%f%f%f%f','headerlines',4);
    
    fid = fopen([folder_name filename]);
    data = textscan(fid,'%f%f%f%f%f%f%f','headerlines',4);
    fclose(fid);
    periph_signal=data{1,1};
    central_signal=data{1,2};
    periph_pulse=data{1,3};
    central_pulse=data{1,4};
    periph_pulse=periph_pulse(~isnan(periph_pulse));     % remove NaNs
    central_pulse=central_pulse(~isnan(central_pulse));     % remove NaNs
end

% % open waitbar
% if record_no==1
%     h = waitbar(record_no/no_of_files,'Processing files...');
% else
%     waitbar(record_no/no_of_files,h,'Processing files...');
% end


% call function to fit reservoir to radial data
[P_av, Pr_av,Pinf_av,Pn_av,Tn_av,Sn_av, fita_av, fitb_av,rsq_av]=fitres_v6(periph_pulse,sampling_rate);
% calculate mean arterial pressure from waveform
%map=str2double(spdata{2,32}); % data in file

% duration of diastole
Tdia=(length(P_av)-Sn_av)/sampling_rate;

% make times for plots
Time=1/sampling_rate*(0:(length(P_av)-1));
TimeR=1/sampling_rate*(0:(length(Pr_av)-1));  % added to allow for different lengths of p_av and Pr_av due to deletion of upturns
%T=1/sampling_rate*(0:(length(periph_pulse)-1));

% Make some calculations using the average data
P_less_dia=P_av-min(P_av);
Pr_less_dia=Pr_av-min(Pr_av);
Pn_less_dia=Pn_av-min(P_av);
Pxs=P_av-Pr_av;

% more filenames
%csvname = regexprep(filename,'.txt','res.csv');    % for excel file - not being used
%csvnamef = fullfile(folder_name,csvname);
id_2 = regexprep(filename,'.txt','');

%% make figures and data subfolders
figfolder='C:\Spdata\figures\';
datafolder='C:\Spdata\results\';
if ~exist(figfolder, 'dir')
    mkdir(figfolder);
end
if ~exist(datafolder, 'dir')
    mkdir(datafolder);
end
% save separate jpg (for viewing) and wmf (for editing)
expression = 'txt';
replace1 = 'wmf';
replace2 ='jpg';
wmffile = regexprep(filename,expression,replace1);
jpgfile = regexprep(filename,expression,replace2);
%%%%%%%%%%%%%

%% perform reservoir analysis on central pressure waveform as a prelude to estimating Pb and Pf
[cP_av,cPr_av,cPinf_av, cPn_av, cTn_av,cSn_av, cfita_av,...
    cfitb_av, crsq_av]=fitres_v6(central_pulse,sampling_rate);

% duration of diastole
cTdia=(length(cP_av)-cSn_av)/sampling_rate;

% calculate central excess pressure
cPxs_av=cP_av-cPr_av;
% Create estimates of Pf and Pb using the assumption that Pb = cPr/2
cPb_av=(cPr_av-min(cP_av))/2;
cPf_av=cP_av-cPb_av-min(cP_av);
Pb_Pf=max(cPb_av)/max(cPf_av);
RI=max(cPb_av)/(max(cPf_av)+ max(cPb_av));

%% wave intensity analysis (currently only done on central P)
% tweak the pressure waveforms cP_av and cPxs_Av to make the wave
% intensity plot look normal (Sphygmocor data starts at the foot
% or sometimes earlier!) - so, there isnt much 'lead in' on the wave intensity if this isnt done.
% But it occasionaly creates a bit of a hiccup at the start.
cpwia = zeros(1,length(cP_av));
cpwia(1:11)=cP_av(end-10:end);
cpwia(12:end)=cP_av(1:end-11);
cpwia=(cpwia*mmHgPa);
cuxwia = zeros(1,length(cPxs_av));
cuxwia(1:11)=cPxs_av(end-10:end);
cuxwia(12:end)=cPxs_av(1:end-11);
% convert xwia to flow velocity and assume peak velocity = 1m/s
% based on Lindroos, M., et al. J Am Coll Cardiol 1993; 21: 1220-1225.
% cuxwia=uconst*(cuxwia-min(cuxwia))/(max(cuxwia)-min(cuxwia));
cuxwia=uconst*(cuxwia/max(cuxwia));

% calculate derivatives
[b,g]=sgolay(Npoly,Frame);   % Calculate S-G coefficients
HalfWin=((Frame+1)/2) -1;
%     p=P-min(P); u=U;         % retain P and U as the original data
N=length(cpwia);
dp=zeros(N,1); duxs=dp;
for n=(Frame+1)/2:N-(Frame+1)/2
    % Zero-th derivative (smoothing only)
    %     ps(n)=dot(g(:,1),p(n-HalfWin:n+HalfWin));
    %     us(n)=dot(g(:,1),u(n-HalfWin:n+HalfWin));
    % 1st differential
    dp(n)=dot(g(:,2),cpwia(n-HalfWin:n+HalfWin));    % pressure difference
    duxs(n)=dot(g(:,2),cuxwia(n-HalfWin:n+HalfWin)); % velocity difference
end

di=dp.*duxs;
di=di*length(dp)^2;             % units fixed - now in W/m2 per cycle^2
minpeak=max(di)/20;             % peaks <5% of Wf1 ignored
% I've left the warning when there are no peaks for now but if it
% needs to be suppressed then 'signal:findpeaks:largeMinPeakHeight' is
% its id
[~,lsys]=min(dp);               % restrict analysis to systole
lsys=round(lsys)+5;             % round and add 5 samples to give a margin for error for duration of systole
[dippks(1),diplocs(1), dipw(1)]=findpeaks(di(1:lsys), 'NPeaks',1,'MinPeakHeight',minpeak); % find 1st dI+ peaks (Wf1)
[dimpks,dimlocs,dimw]=findpeaks(-di(1:lsys), 'NPeaks',1,'MinPeakHeight',0.7*max(-di)); % find one dI- peaks
% to find Wf2 flip the di in systole and find peak.
[dippks(2),diplocs(2), dipw(2)]=findpeaks(flipud(di(1:lsys)), 'NPeaks',1,'MinPeakHeight',minpeak); % find 2nd dI+ peaks (Wf2) by flipping the data and running findpeak in the other direction
diplocs(2)=lsys-diplocs(2)+1;     % correct location for the flipping
% check peaks
%     figure; hold on; plot(di); plot(diplocs(1),dippks(1),'ko');
%     plot(dimlocs,-dimpks,'ro'); plot(diplocs(2),dippks(2),'ks');
dipt=diplocs/sampling_rate;
dimt=dimlocs/sampling_rate;

% error trap for when Wf2 is unmeasureable - with new routine this may
% be unecessary
if length(dippks)==1
    dippks(2)=0;
    dipt(2)=0;
    diparea(2)=0;
end

% calculate areas
% For a Gaussian curve (assumed) the area is 1.06447*height*width
diparea=1.06447*dippks.*dipw;
dimarea=1.06447*dimpks.*dimw;
wri=dimarea/diparea(1);

% Estimate c (wavespeed) as k*dP/du where k is empirical constant
% currently k = 1!
rhoc=max(Pxs)*mmHgPa/1000; % fixed units (m/s)

%% HRV
%% perform HRV analysis on peripehral pressure waveform
[P_all, sdsbp, nbeats, rr_ms,rrS_ms, sdnn_ms, ...
    sdnnS_ms, rmssd_ms, rmssdS_ms, brs_ms_mmhg,sysloc,dialoc]...
    =hrv_v1(periph_pulse,central_signal,sampling_rate);

% perform HRV analysis on central pressure waveform
[cP_all, csdsbp, cnbeats, crr_ms,crrS_ms, csdnn_ms, ...
    csdnnS_ms, crmssd_ms, crmssdS_ms, cbrs_ms_mmhg,csysloc, cdialoc]...
    =hrv_v1(central_pulse,central_signal,sampling_rate);
%% Print figures
% Pressure and central reservoir
T_all=1/sampling_rate*(0:(length(P_all)-1));
TimeC=(1:length(cP_av))/sampling_rate;
%clf
% figure('visible','off');                     % dont display figure
figure('visible','on');  
subplot(1,2,1); hold on;
plot(T_all,P_all);plot(sysloc/sampling_rate, P_all(sysloc),'ro');
plot(dialoc/sampling_rate, P_all(dialoc),'rs');
xlabel('Time (s)')
ylabel('BP (mmHg)')
title('Pulse traces')
box off;
subplot(1,2,2);
plot(TimeC,cP_av,TimeC,cPr_av,'r', TimeC,cPxs_av,'k');
xlabel('Time (s)')
ylabel('BP (mmHg)')
title('P, Pres, Pxs')
box off;
% % print ('-dmeta', '-r300' , [figfolder wmffile]);
% print ('-djpeg', '-r300' , [figfolder jpgfile]);

% WI
%clf
% figure('visible','off');                     % dont display figure
figure('visible','on');  
TimeDI=(1:length(di))/sampling_rate;
plot (TimeDI, di);  hold on;                             % **** to allow for new length
plot(dipt(1),dippks(1),'ko');
plot(dimt,-dimpks,'ro'); plot(dipt(2),dippks(2),'ks');
xlabel('Time (s)')
ylabel('dI (W/m^2/cycle^2)')
title('Wave intensity with Wf1, Wb, Wf2 identified')
box off;
% wmffile1 = regexprep(filename,'.txt','w.wmf');
% jpgfile1 = regexprep(filename,'.txt','w.jpg');
% % print ('-dmeta', '-r300' , [figfolder wmffile1]);
% print ('-djpeg', '-r300' , [figfolder jpgfile1]);
drawnow();                  % added to attempt to stop java leak

% P, Pf, Pb
%clf
% figure('visible','off');                     % dont display figure
figure('visible','on');  
hold on; plot(TimeC, cPf_av,'b', TimeC, cPb_av,'r', TimeC, cP_av-min(cP_av),'k');
xlabel('Time (s)')
ylabel('BP (mmHg)')
title('P, Pf and Pb')
% wmffile2 = regexprep(filename,'.txt','fb.wmf');
% jpgfile2 = regexprep(filename,'.txt','fb.jpg');
% % print ('-dmeta', '-r300' , [figfolder wmffile2]);
% print ('-djpeg', '-r300' , [figfolder jpgfile2]);

% message at end
msg=string(no_of_files);
msg=strcat(msg, " file(s) processed");
msgbox(msg, 'Done');
