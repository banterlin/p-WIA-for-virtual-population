%% Script adapted from main_v06 to fit Stroke Volumn impaired patients
% 
%% Raw dataset can be obtained from https://github.com/ryanreavette/WaveIntensityData
%% batch_res_v14 - batch analysis of radial pulse data using kreservoir_v14
%% Copyright 2020 Alun Hughes & Kim Parker
%% Modified by Lin for Project evluating effectiveness of pressure-only WIA against normal WIA.
% This software is distributed under under the terms of the GNU General Public License
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

%%%%%%%%%%%%%%%
%% m files required to be in directory
% fitres_v6.m
% kreservoir_v14.m
% sys_dia_finder.m
% WIA.m
% ReadFileNames.m
%%%%%%%%%%%%%%%%

%% A differnt data aqusition method is used here, where folder directory is obtained first and then we analyse
%% Constants
clear all; clc;

close all;


sampling_rate = 1000;  % For our sample interval of 1ms.

kres_v=14;          % Version tracking for reservoir fitting
headernumber=28;       % headers for columns of results (see end)
mmHgPa = 133;          % P conversion for WIA
inlet = 0;              % If we are dealing with inlet data, inlet=1. (0 otherwise)
uconst=1;              % empirical constant to convert normalized
% velocity to m/s based on Hughes et al.
% Front. Physiol. DOI: 10.3389/fphys.2020.00550
Npoly=3;               % Order of polynomial fit for sgolay
Frame=61;               % Window length for sgolay
version='0.6';         % Version of Main
%%%%%%%%%%%%%%%%

[ FList ] = ReadFileNames('C:\Spdata\SV_impaired\population\controls\');
Index = strfind(FList, 'right-common-carotid');
Idx = find(~cellfun(@isempty,Index));
%% Loop through all dataset

proc_var=cell(length(Idx),headernumber);
record_no=1;
for file_number = 1:length(Idx)
    Filename = FList(Idx(file_number));
    record_no = file_number;
    
    fid = fopen(char(Filename));
    data = textscan(fid,'%f%f%f%f%f','headerlines',1);
    fclose(fid);
        pressure=data{1,4}/mmHgPa*100; % RAW Pressure arranged in Pa/100, converting this to mmHg 
    %     central_pulse=data{1,4};
    velocity = data{1,3}./100;  % Convert raw data cm/s to m/s
% Following code checks the peaks between 1 detected cycle with 1.5
% detected cycle, to avoid cases where young subjects showed very low
% magnitude of dirotic notch (where the algorithm would see it as a diastolic trough)
    offset = 3500;
    [sysloc, dialoc]=sys_dia_finder(pressure(offset:end),sampling_rate) ;
    dialoc = dialoc + offset;
    if max(pressure(dialoc(end-2) : dialoc(end-1))) - max(pressure(dialoc(end-3) : dialoc(end-2))) >2 ...
            && dialoc(end-1) - dialoc(end-2) < 800
        periph_T  = [dialoc(end-4) : dialoc(end-2)];
    elseif max(pressure(dialoc(end-2) : dialoc(end-1))) - max(pressure(dialoc(end-3) : dialoc(end-2))) < -2 ...
            && dialoc(end-1) - dialoc(end-2) < 800
        periph_T  = [dialoc(end-3) : dialoc(end-1)];
    else
        periph_T  = [dialoc(end-2) : dialoc(end-1)];
    end

    
    single_pulse = pressure(periph_T,1);    % Single cardiac cycle analysis
    velocity = velocity(periph_T,1);        % Single cardiac cycle analysis
    
    single_pulse=single_pulse(~isnan(single_pulse));     % remove NaNs
%     central_signal= central_signal(~isnan(central_signal));
    velocity = velocity (~isnan(velocity));
    %     central_pulse=central_pulse(~isnan(central_pulse));     % remove NaNs
    
    % call function for reservoir fitting
    [P_av, Pr_av,Pinf_av,Pn_av,Tn_av,Sn_av, fita_av, fitb_av,rsq_av,dp,nn]=fitres_v6(single_pulse,sampling_rate,kres_v);
    
    Pxs=P_av-Pr_av; % Excess Pressure
    
    zci = @(v) find(diff(sign(v))); % Function for finding zero crossing
    n_Zx = length(zci(dp));         % Number of zero crossing presented by the pressure wave
    %% Diastolic Excess pressure zeroing.  -- designed for inlet processing.
    if inlet == 1;
        idx1= zci(Pxs);
        if length(idx1) > 2 ;
            idx2 = idx1(3);
        else
            idx2 = idx1(2);
        end
        Pxs(idx2:end) = 0;
    end
%     switch answer
%         case 'Yes'
%             folder_name = uigetdir;
%         case 'No [end]'             % end if no folder identified
%             return
%     end
% end

    %% Print Figure for reservior analysis
    
    TimeC=(1:length(P_av))/sampling_rate;

%     figure('visible','off');                     % dont display figure
    

%     f1 = figure(1); set(f1,'Color','w');
%     plot(TimeC,P_av,TimeC,Pr_av,'r', TimeC,Pxs,'k');legend('P','Pres','Pxs')
%     xlabel('Time (s)')
%     ylabel('BP (mmHg)')
%     title('P, Pres, Pxs')
%     box off;
    
    % calculate mean arterial pressure from waveform
    %map=str2double(spdata{2,32}); % data in file
    % calculate mean arterial pressure from waveform
    cPb_av=(Pr_av-min(P_av))/2;     % Half of Pres elevation
    % Assuming backward traveling wave is half of the reservior pressure
    cPf_av=P_av-cPb_av-min(P_av);
    Pb_Pf=max(cPb_av)/max(cPf_av);
    RI=max(cPb_av)/(max(cPf_av)+ max(cPb_av));
    %% Wave intensity analysis (WI Using pressure-only method)
    % tweak the pressure waveforms cP_av and cPxs_Av to make the wave
    % intensity plot look normal (Sphygmocor data starts at the foot
    % or sometimes earlier!) - so, there isnt much 'lead in' on the wave intensity if this isnt done.
    % But it occasionaly creates a bit of a hiccup at the start.
    cpwia = zeros(1,length(P_av));
    P_wi = P_av;
    cpwia(1:11)=P_wi(end-10:end);
    cpwia(12:end)=P_wi(1:end-11);
    % cpwia = P_av;

    cuxwia = zeros(1,length(Pxs));
    cuxwia(1:11)=Pxs(end-10:end);
    cuxwia(12:end)=Pxs(1:end-11);
    % convert xwia to flow velocity and assume peak velocity = 1m/s
    % based on Lindroos, M., et al. J Am Coll Cardiol 1993; 21: 1220-1225.
    % cuxwia=uconst*(cuxwia-min(cuxwia))/(max(cuxwia)-min(cuxwia));
    
%     uconst = max(velocity); % Using actual max(U) for better estimation
    cuxwia=uconst*(cuxwia/max(cuxwia));   % Peak Velocity calibrated to 1m/s assumption/max(U) depend on uconst


    % Estimate c (wavespeed) as k*dP/du where k is empirical constant
    % currently k = 1!
    rhocxs=max(Pxs)./1060./uconst; % fixed units (kg) 
    [dixs,wrixs,diplocsxs,dippksxs,dimlocsxs,dimpksxs,diptxs,dimtxs,dipareaxs,dimareaxs] = WIA(Npoly,Frame,cpwia,cuxwia,sampling_rate);
   
    %% Wave intensity analysis(True WI calculation)
    P = single_pulse;
    p=P-min(P);
    p = zeros(1,length(P_av));
    p(1:11)=P_av(end-10:end);
    p(12:end)=P_av(1:end-11);

    u = zeros(1,length(velocity));
    u(1:11)=velocity(end-10:end);
    u(12:end)=velocity(1:end-11);
    % u=uconst*(u/max(u));    %Make velocity unitless
    % u=u/100;       % retain P and U as the original data
    wavespeed = max(data{1,5}./100);
    [di,wri,diplocs,dippks,dimlocs,dimpks,dipt,dimt,diparea,dimarea] = WIA(Npoly,Frame,p,u,sampling_rate);
    
    
    %% Corrcoef between Pxs and velocity from raw data.
    corr_Pxs = corrcoef(velocity(1:length(Pxs)),Pxs);
    % fprintf("Correlation coefficent between Pxs and Raw Velocity is %f \n",corr_Pxs(1,2));
    corrdI = corrcoef(di,dixs);
%     fprintf("Correlation coefficent between dI and P-only dI is %f \n",corrdI(1,2));
%     fprintf("r^2 of diastolic fit is %f \n",rsq_av);
    % Comparing wave intensities
    % Constructing time interval over sampling period
    t = [0:length(di)-1]/sampling_rate;
    

%     f2 = figure(2);set(f2,'Color','w');
%     subplot(1,2,1);plot(t,di);title('P-U WI');
%     xlabel('time(m/s)');ylabel('Wave intensity(W/m^2)')
%     subplot(1,2,2);plot(t,dixs);title('Pressure-only WI');
%     xlabel('time(m/s)');ylabel('Wave intensity(W/m^2)')
%     txt = ['Correlation:',num2str(corrdI(1,2))];sgtitle(txt);
    % Comparing velocity with Pxs
    t = [0:length(velocity)-1]/sampling_rate;
    t2 = [0:length(Pxs)-1]/sampling_rate;

%     f3 = figure(3);set(f3,'color','w'); 
%     subplot(1,2,1);plot(t,velocity);title('True Velocity waveform');
%     xlabel('time(m/s)');ylabel('Velocity(m/s)');
%     subplot(1,2,2);plot(t2,Pxs);title('Pxs');
%     xlabel('time(m/s)');ylabel('Velocity(m/s)');
%     txt = ['Correlation:',num2str(corr_Pxs(1,2))];sgtitle(txt);
    
    %% Variable calculations

    [max_p_av,max_tp_av]=max(P_av);
    [max_pxs,max_txs]=max(Pxs); % calculation
    [max_pr,max_tr]=max(Pr_av); % calculation
    Pr_less_dia=Pr_av-min(Pr_av);
    [max_prld,max_trld]=max(Pr_less_dia); % calculation
    
    SBP = max(P_av);
    DBP = min(P_av);
    AP = (SBP - P_av(nn));
    AIx = ((SBP - P_av(nn)) / (SBP-DBP))*100;

%     Pr_less_dia=Pr_av-min(Pr_av);
%     [max_prld,max_trld]=max(Pr_less_dia); % calculation
    euc_dist = dtw(di,dixs);        % DTW distance 
    max_velocity = max(velocity);
    %% write variables
    proc_var{record_no,1}= Filename;         % filename
    proc_var{record_no,2}= dippks(1);       % Wf1 intensity
    proc_var{record_no,3}= dimpks;          % Wb intensity
    proc_var{record_no,4}= dippks(2);       % Wf2 intensity
    proc_var{record_no,5}= dippksxs(1);         % Timing of Wf1 wave
    proc_var{record_no,6}= dimpksxs;         % Timing of Wf1xs wave
    proc_var{record_no,7}= dippksxs(2);          % Conversion used as constant velocity assumption
    proc_var{record_no,8}= dimlocs;          % Conversion used as constant velocity assumption
    proc_var{record_no,9}= dimlocsxs;          % Conversion used as constant velocity assumption
    proc_var{record_no,10}=sampling_rate;   % sampling rate, 1/sec [used to be time of max Pres-diastolic but since this == Time max Pres I replaced it with the sampling rate]
    proc_var{record_no,11}=max_pxs;         % max excess P, mmHg
    proc_var{record_no,12}=Tn_av;           % Time of end systole by max -dp/dt, sec
    proc_var{record_no,13}=Pinf_av;         % P infinity, mmHg
    proc_var{record_no,14}=Pn_av;           % P at end systole by max -dp/dt, mmHg
    proc_var{record_no,15}=wri;         % True Wave Reflection index
    proc_var{record_no,16}=wrixs;         % Excess Wave Relfection Index
    proc_var{record_no,17}=Pb_Pf;           % Pb/Pf (ratio of backward to forward pressure - also know as reflection magnitude. RM
    proc_var{record_no,18}=diplocs(1);
    proc_var{record_no,19}=diplocs(2);
    proc_var{record_no,20}=diplocsxs(1);
    proc_var{record_no,21}=diplocsxs(2);
    proc_var{record_no,22}=corr_Pxs(1,2);          % PCC of Pxs with velocity
    proc_var{record_no,23}=corrdI(1,2);          % PCC of dIxs with dI
    proc_var{record_no,24}=euc_dist;          % DTW distance
    proc_var{record_no,25}=rhocxs;              % Estimated Wavespeed
    proc_var{record_no,26}=max(data{1,5})/100;  % Actual Wavespeed
    proc_var{record_no,27}=max(velocity)/100;  % Maximun velocity
    proc_var{record_no,28}=version;            % Version of reservior analysis used

    

    fprintf('Processing............... %f%% \n', file_number/length(Idx).*100)
end
%% 
datafolder='C:\Spdata\results\';
if ~exist(datafolder, 'dir')
    mkdir(datafolder);
end
xlsfile='C:\Spdata\results\resdata.xlsx'; % Or use .xls for compatibility

header = {'patid' 'dippks1' 'dimpks' 'dippks2' 'dippksxs1' 'dimpksxs'...
    'dippksxs2','dimlocs','dimlocsxs',	're_sam_rate', 're_maxxsp'...
    're_tn' 're_pinf' 're_pn' 'wri' 'wrixs','re_pb_pf', ...
    'diploc1','diploc2','diplocxs1','diplocxs2','corr_Pxs', 'corr_dI', ... 
    'DTW_dist','rhocxs','wavespeed','U_max''kres_v'}; % header
% header = {'re_file' 'dippks1' 'dimpks' 'dippks2' 'dippksxs1' 'dimpksxs'...
%     'dippksxs2','dimlocs','dimlocsxs',	're_sam_rate', 're_maxxsp'...
%     're_tn' 're_pinf' 're_pn' 'wri' 'wrixs','re_pb_pf', ...
%     'diploc1','diploc2','diplocxs1','diplocxs2','corr_Pxs', 'corr_dI', ... 
%     'DTW_dist','rhocxs','wavespeed','dI+_corr','dI-_corr','kres_v','n_Zx'}; % header

% writetable
Results_table=cell2table(proc_var, 'VariableNames',header);
writetable(Results_table, xlsfile);

clc;
fprintf('Analysis Complete')


