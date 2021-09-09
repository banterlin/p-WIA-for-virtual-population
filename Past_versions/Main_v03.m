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

%% Versions modified by Lin
%  v0.1 (First version working for virtual population with pulse and signal input)
%  v0.2 (** Removed code for calculating variables that are not considered in our study
%           Wave intensity calculation from Pxs.
%           automatic detection of last cardiac cycle from input of virtual population
%  v0.3  Separated WI calculation into the function WIA. added results
%  export to excel
%%%%%%%%%%%%%%%
%% m files required to be in directory
% fitres_v6.m
% kreservoir_v14.m
% hrv_v1.m
%%%%%%%%%%%%%%%%
%% Constants
clear all; clc;
% close all;

sampling_rate = 1000;  % For our sample interval of 1ms.

kres_v='v14';          % Version tracking for reservoir fitting
headernumber=22;       % headers for columns of results (see end)
mmHgPa = 133;          % P conversion for WIA
uconst=1;              % empirical constant to convert normalized
% velocity to m/s based on Hughes et al.
% Front. Physiol. DOI: 10.3389/fphys.2020.00550
Npoly=3;               % Order of polynomial fit for sgolay
Frame=61;               % Window length for sgolay
version='0.3';         % Version of Main
%%%%%%%%%%%%%%%%
%% Select files
% default folder as per manual
folder_name ='C:\Spdata\single\';
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
    record_no = file_number;        %for multiple file input
    filename=file_lists(record_no).name;
    
    fid = fopen([folder_name filename]);
    data = textscan(fid,'%f%f%f%f%f','headerlines',5);
    fclose(fid);
    %     periph_signal=data{1,1};
%     central_signal=data{1,2};
    single_pulse=data{1,3};
    %     central_pulse=data{1,4};
    velocity = data{1,1}./100;  % Convert raw data cm/s to m/s
    
    single_pulse=single_pulse(~isnan(single_pulse));     % remove NaNs
%     central_signal= central_signal(~isnan(central_signal));
    velocity = velocity (~isnan(velocity));
    %     central_pulse=central_pulse(~isnan(central_pulse));     % remove NaNs
    
    
    % call function for reservoir fitting
    [P_av, Pr_av,Pinf_av,Pn_av,Tn_av,Sn_av, fita_av, fitb_av,rsq_av]=fitres_v6(single_pulse,sampling_rate);
    
    % calculate mean arterial pressure from waveform
    %map=str2double(spdata{2,32}); % data in file
    cPb_av=(Pr_av-min(P_av))/2;
    cPf_av=P_av-cPb_av-min(P_av);
    Pb_Pf=max(cPb_av)/max(cPf_av);
    RI=max(cPb_av)/(max(cPf_av)+ max(cPb_av));
    
    % duration of diastole
    % Tdia=(length(P_av)-Sn_av)/sampling_rate;
    
    Pxs=P_av-Pr_av; % Excess Pressure
    
    
    %% Wave intensity analysis (WIUsing pressure-only method)
    % tweak the pressure waveforms cP_av and cPxs_Av to make the wave
    % intensity plot look normal (Sphygmocor data starts at the foot
    % or sometimes earlier!) - so, there isnt much 'lead in' on the wave intensity if this isnt done.
    % But it occasionaly creates a bit of a hiccup at the start.
    cpwia = zeros(1,length(P_av));
    cpwia(1:11)=P_av(end-10:end);
    cpwia(12:end)=P_av(1:end-11);
    % cpwia = P_av;
    cpwia=(cpwia*mmHgPa);
    cuxwia = zeros(1,length(Pxs));
    cuxwia(1:11)=Pxs(end-10:end);
    cuxwia(12:end)=Pxs(1:end-11);
    % convert xwia to flow velocity and assume peak velocity = 1m/s
    % based on Lindroos, M., et al. J Am Coll Cardiol 1993; 21: 1220-1225.
    % cuxwia=uconst*(cuxwia-min(cuxwia))/(max(cuxwia)-min(cuxwia));
    cuxwia=uconst*(cuxwia/max(cuxwia));
    
    [dixs,wrixs,diplocsxs,dippksxs,dimlocsxs,dimpksxs] = WIA(Npoly,Frame,cpwia,cuxwia,sampling_rate);
    
    
    % Estimate c (wavespeed) as k*dP/du where k is empirical constant
    % currently k = 1!
%     rhocxs=max(Pxs)*mmHgPa/1000; % fixed units (m/s)
    
    %% Wave intensity analysis(True WI calculation)
    P = single_pulse;
    p=P-min(P);
    p = zeros(1,length(P_av));
    p(1:11)=P_av(end-10:end);
    p(12:end)=P_av(1:end-11);
    p=(p*mmHgPa);
    u = zeros(1,length(velocity));
    u(1:11)=velocity(end-10:end);
    u(12:end)=velocity(1:end-11);
    % u=uconst*(u/max(u));    %Make velocity unitless
    % u=u/100;       % retain P and U as the original data
    
    [di,wri,diplocs,dippks,dimlocs,dimpks] = WIA(Npoly,Frame,p,u,sampling_rate);
    
    %% Corrcoef between Pxs and velocity from raw data.
    corr_Pxs = corrcoef(velocity(1:length(Pxs)),Pxs);
    % fprintf("Correlation coefficent between Pxs and Raw Velocity is %f \n",corr_Pxs(1,2));
    corrdI = corrcoef(di,dixs);
    % fprintf("Correlation coefficent between dI and P-only dI is %f \n",corrdI(1,2));
    % fprintf("r^2 of diastolic fit is %f \n",rsq_av);
    %% Comparing wave intensities
    % Constructing time interval over sampling period
    t = [0:length(di)-1]/sampling_rate;
    f2 = figure(2);set(f2,'Color','w'); subplot(1,2,1);plot(t,di);title('Wave intensity from conventional method');
    xlabel('time(m/s)');ylabel('Wave intensity(W/m^2)')
    subplot(1,2,2);plot(t,dixs);title('Wave intensity from pressure-only method');
    xlabel('time(m/s)');ylabel('Wave intensity(W/m^2 *10^4 / Cycle^2)')
    %% Comparing velocity with Pxs
    t = [0:length(velocity)-1]/sampling_rate;
    t2 = [0:length(Pxs)-1]/sampling_rate;
    f82 = figure(82);set(f82,'Color','w');
    subplot(1,2,1);plot(t,velocity);title('True Velocity waveform');
    xlabel('time(m/s)');ylabel('Velocity(m/s)');
    subplot(1,2,2);plot(t2,Pxs);title('Pxs');
    xlabel('time(m/s)');ylabel('Velocity(m/s)');
    
    
    
    %% HRV
    %% perform HRV analysis on peripehral pressure waveform
    % [P_all, sdsbp, nbeats, rr_ms,rrS_ms, sdnn_ms, ...
    %     sdnnS_ms, rmssd_ms, rmssdS_ms, brs_ms_mmhg,sysloc,dialoc]...
    %     =hrv_v1(periph_pulse,central_signal,sampling_rate);
    
    
%     [P_all,sdsbp ,sysloc, dialoc]=hrv(single_pulse,central_signal,sampling_rate);
    
    
    %% Print figures
    % Pressure and central reservoir
%     T_all=1/sampling_rate*(0:(length(P_all)-1));
    % TimeC=(1:length(cP_av))/sampling_rate;
    TimeC=(1:length(P_av))/sampling_rate;
    %clf
    % figure('visible','off');                     % dont display figure
    f1 = figure(1); set(f1,'Color','w');
    % subplot(1,2,1); hold on;
    % plot(T_all,P_all);plot(sysloc/sampling_rate, P_all(sysloc),'ro');
    % plot(dialoc/sampling_rate, P_all(dialoc),'rs');
    % xlabel('Time (s)')
    % ylabel('BP (mmHg)')
    % title('Pulse traces')
    % box off;
    % subplot(1,2,2);
    plot(TimeC,P_av,TimeC,Pr_av,'r', TimeC,Pxs,'k');
    xlabel('Time (s)')
    ylabel('BP (mmHg)')
    title('P, Pres, Pxs for M60-69')
    box off;
    % % print ('-dmeta', '-r300' , [figfolder wmffile]);
    % print ('-djpeg', '-r300' , [figfolder jpgfile]);
    
    % WI
    %clf
    % figure('visible','off');                     % dont display figure
    % figure('visible','on');
    % TimeDI=(1:length(di))/sampling_rate;
    % plot (TimeDI, di);  hold on;                             % **** to allow for new length
    % plot(dipt(1),dippks(1),'ko');
    % plot(dimt,-dimpks,'ro'); plot(dipt(2),dippks(2),'ks');
    % xlabel('Time (s)')
    % ylabel('dI (W/m^2/cycle^2)')
    % title('Wave intensity with Wf1, Wb, Wf2 identified')
    % box off;
    % wmffile1 = regexprep(filename,'.txt','w.wmf');
    % jpgfile1 = regexprep(filename,'.txt','w.jpg');
    % print ('-dmeta', '-r300' , [figfolder wmffile1]);
    % print ('-djpeg', '-r300' , [figfolder jpgfile1]);
    % drawnow();                  % added to attempt to stop java leak
    
    % P, Pf, Pb
    %clf
    % figure('visible','off');                     % dont display figure
    % figure('visible','on');
    % hold on; plot(TimeC, cPf_av,'b', TimeC, cPb_av,'r', TimeC, cP_av-min(cP_av),'k');
    % xlabel('Time (s)')
    % ylabel('BP (mmHg)')
    % title('P, Pf and Pb')
    % wmffile2 = regexprep(filename,'.txt','fb.wmf');
    % jpgfile2 = regexprep(filename,'.txt','fb.jpg');
    % % print ('-dmeta', '-r300' , [figfolder wmffile2]);
    % print ('-djpeg', '-r300' , [figfolder jpgfile2]);
    
    
    % %% message at end
    % msg=string(no_of_files);
    % msg=strcat(msg, " file(s) processed");
    % msgbox(msg, 'Done');
    
    
    %% Variable calculations

    [max_p_av,max_tp_av]=max(P_av);
    [max_pxs,max_txs]=max(Pxs); % calculation
    [max_pr,max_tr]=max(Pr_av); % calculation
%     Pr_less_dia=Pr_av-min(Pr_av);
%     [max_prld,max_trld]=max(Pr_less_dia); % calculation

    
    
    %% write variables
    proc_var{record_no,1}=filename;         % filename
    proc_var{record_no,2}=max_p_av;         % max P (SBP), mmHg
    proc_var{record_no,3}=max_tp_av/sampling_rate; % time max P (SBP), s
    proc_var{record_no,4}=min(P_av);        % min P, DBP, mmHg
    proc_var{record_no,5}=sum(Pr_av)/sampling_rate; % integral Pres, mmHg.s
    proc_var{record_no,6}=max_pr;           % max Pres, mmHg
    proc_var{record_no,7}=max_tr/sampling_rate; % Time to max Pres, sec
%     proc_var{record_no,8}=sum(Pr_less_dia)/sampling_rate; % integral Pres-diastolic, mmHg.s
%     proc_var{record_no,9}=max_prld;         % max Pres-diastolic, mmHg
    proc_var{record_no,8}=sampling_rate;   % sampling rate, 1/sec [used to be time of max Pres-diastolic but since this == Time max Pres I replaced it with the sampling rate]
    proc_var{record_no,9}=sum(Pxs)/sampling_rate; % Integral excess pressure, mmHg.s
    proc_var{record_no,10}=max_pxs;         % max excess P, mmHg
    proc_var{record_no,11}=max_txs/sampling_rate; % time of max excess P, sec
    proc_var{record_no,12}=Tn_av;           % Time of end systole by max -dp/dt, sec
    proc_var{record_no,13}=Pinf_av;         % P infinity, mmHg
    proc_var{record_no,14}=Pn_av;           % P at end systole by max -dp/dt, mmHg
    proc_var{record_no,15}=fita_av;         % rate constant A (ka), 1/sec
    proc_var{record_no,16}=fitb_av;         % rate constant B (kb), 1/sec
    proc_var{record_no,17}=rsq_av;          % R2 for diastolic fit
    proc_var{record_no,18}=Pb_Pf;           % Pb/Pf (ratio of backward to forward pressure - also know as reflection magnitude. RM
    proc_var{record_no,19}=RI;              % Reflection index (RI) calculated as Pb/(Pb+Pf)
%     proc_var{record_no,20}=dippks(1);       % W1 intensity
%     proc_var{record_no,36}=dipt(1);         % W1 time
%     proc_var{record_no,37}=diparea(1);      % W1 area
%     proc_var{record_no,21}=dimpks;          % W-1 intensity
%     proc_var{record_no,39}=dimt;            % W-1 time
%     proc_var{record_no,40}=dimarea;         % W-1 area
%     proc_var{record_no,22}=dippks(2);       % W2 intensity
%     proc_var{record_no,42}=dipt(2);         % W2 time
%     proc_var{record_no,43}=diparea(2);      % W2 area
%     proc_var{record_no,44}=wri;             %
    proc_var{record_no,20}=wrixs;            % WRI for excess pressure
%     proc_var{record_no,46}=rhocxs;              % rhoc estimation for reservoir analysis
    proc_var{record_no,21}=kres_v;          % version tracking
    proc_var{record_no,22}=version;         % version of main used
    
    
end

%% Save the results as an excel spreadheet
% datafolder='C:\Spdata\results\';
% if ~exist(datafolder, 'dir')
%     mkdir(datafolder);
% end
% xlsfile='C:\Spdata\results\resdata.xls';
% header = {'re_file' 're_maxp' 're_tmaxp' 're_minp'	're_intpr' 're_maxpr'...
%     're_tmaxpr'	're_sam_rate', 're_intxsp','re_maxxsp'...
%     're_tmaxxsp' 're_tn' 're_pinf' 're_pn' 're_fita' 're_fitb'	're_rsq' ...
%     're_pb_pf', 're_ri', 're_wri','re_version','version'}; % header
% 
% % writetable
% Results_table=cell2table(proc_var, 'VariableNames',header);
% writetable(Results_table, xlsfile);


