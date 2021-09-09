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
%  v0.4 added batch analysis to all the txt files in Spata folder.
%  v0.5 optimised single pressure waveform detection --> imporved results
%       especially for the younger popoulations
%  v05 is the latest working version as of 04/31
%  V06 enables automatic variable expoential fitting time for the reservoir
%       pressure analysis
%%%%%%%%%%%%%%%
%% m files required to be in directory
% fitres_v6.m
% kreservoir_v14.m
% sys_dia_finder.m
% WIA.m
%%%%%%%%%%%%%%%%

clear all; clc;
% close all;
%% Constants
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
%% Select files
% default folder as per manual
folder_name ='C:\Spdata\LCC\';
% check that C:\Spdata exists and if not allows new folder to be chosen
if ~exist('C:\Spdata\LCC\', 'dir')
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
%     pressure=data{1,4};
    pressure=data{1,4}/mmHgPa*100; % RAW Pressure arranged in Pa/100, converting this to mmHg 

    %     central_pulse=data{1,4};
    velocity = data{1,3}./100;  % Convert raw data cm/s to m/s
% Following code checks the peaks between 1 detected cycle with 1.5
% detected cycle, to avoid cases where young subjects showed very low
% magnitude of dirotic notch (where the algorithm would see it as a diastolic trough)
    [sysloc, dialoc]=sys_dia_finder(pressure,sampling_rate) ;
    if max(pressure(dialoc(end-2) : dialoc(end-1))) - max(pressure(dialoc(end-3) : dialoc(end-2))) >10
        periph_T  = [dialoc(end-4) : dialoc(end-2)];
    elseif max(pressure(dialoc(end-2) : dialoc(end-1))) - max(pressure(dialoc(end-3) : dialoc(end-2))) < -10
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
    version = kres_v;
    [P_av, Pr_av,Pinf_av,Pn_av,Tn_av,Sn_av, fita_av, fitb_av,rsq_av,dp,nn]=fitres_v6(single_pulse,sampling_rate,version);
    
    Pxs=P_av-Pr_av; % Excess Pressure
    

    
    
%     if max(Pxs(1:400))- max(Pxs(400:end)) < 5; % Threshold Value here require more study
%         version = 16;
%         fprintf('1')
%         [P_av, Pr_av,Pinf_av,Pn_av,Tn_av,Sn_av, fita_av, fitb_av,rsq_av,dp,nn]=fitres_v6(single_pulse,sampling_rate,version);
%         Pxs=P_av-Pr_av; % Excess Pressure
%     end
    

%     f99 = figure(99); set(f99,'Color','w');plot(dp);title('dp')% Figure for rate of change of pressure
    

    zci = @(v) find(diff(sign(v))); % Function for finding zero crossing
    n_Zx = length(zci(dp));         % Number of zero crossing presented by the pressure wave
    %% Diastolic Excess pressure zeroing. 
    if inlet == 1;
        idx1= zci(Pxs);
        if length(idx1) > 2 ;
            idx2 = idx1(3);
        else
            idx2 = idx1(2);
        end
        Pxs(idx2:end) = 0;
    end
    
    
     f2 = figure(36);  set(f2,'Color','w');
     plot(Pxs);
     
    % calculate mean arterial pressure from waveform
    %map=str2double(spdata{2,32}); % data in file
    cPb_av=(Pr_av-min(P_av))/2;
    cPf_av=P_av-cPb_av-min(P_av);
    Pb_Pf=max(cPb_av)/max(cPf_av);
    RI=max(cPb_av)/(max(cPf_av)+ max(cPb_av));
    
    % duration of diastole
    % Tdia=(length(P_av)-Sn_av)/sampling_rate;
    

    
    
    %% Wave intensity analysis (WI Using pressure-only method)
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
    
%     uconst = max(velocity); % Using actual max(U) for better estimation
    cuxwia=uconst*(cuxwia/max(cuxwia));   % Peak Velocity calibrated to 1m/s assumption/max(U) depend on uconst

    
    rhocxs=max(Pxs)*mmHgPa./1060./uconst; % despite the variable name rhoc, this is actually the Estimated Wavespeed
    [dixs,wrixs,diplocsxs,dippksxs,dimlocsxs,dimpksxs,diptxs,dimtxs,dipareaxs,dimareaxs] = WIA(Npoly,Frame,cpwia,cuxwia,sampling_rate);
       
    %% Wave separation Code in this section can be replaced by using WIA_v2.m, instead of WIA.m
%     [~,lsys]=min(dp);               % restrict analysis to systole
%     lsys = lsys + 40;
%     minpeak=max(dixs_plus)/20;             % peaks <5% of Wf1 ignored
%     
%     [dixs_plus_peak(1),dixs_plus_loc(1)]=findpeaks(dixs_plus(1:lsys), 'NPeaks',1,'MinPeakHeight',minpeak); % find 1st dI+ peaks (Wf1)
%     [dixs_plus_peak(2),dixs_plus_loc(2)]=findpeaks(flip(dixs_plus(1:lsys)), 'NPeaks',1,'MinPeakHeight',minpeak); % find 2nd dI+ peaks (Wf2) by flipping the data and running findpeak in the other direction
%     [dixs_minus_peak,dixs_minus_loc]=findpeaks(-dixs_minus(1:lsys), 'NPeaks',1,'MinPeakHeight',0.7*max(-dixs_minus)); % find one dI- peaks
    
    
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
    
    wavespeed = max(data{1,5}./100);
    [di,wri,diplocs,dippks,dimlocs,dimpks,dipt,dimt,diparea,dimarea] = WIA(Npoly,Frame,p,u,sampling_rate);
    
    %% Corrcoef between Pxs and velocity from raw data.
    corr_Pxs = corrcoef(velocity(1:length(Pxs)),Pxs);
    corrdI = corrcoef(di,dixs);
    %% Comparing wave intensities
    % % Constructing time interval over sampling period
    % t = [0:length(di)-1]/sampling_rate;
    % f2 = figure(2);set(f2,'Color','w'); subplot(1,2,1);plot(t,di);title('Wave intensity from conventional method');
    % xlabel('time(s)');ylabel('Wave intensity(W/m^2)')
    % subplot(1,2,2);plot(t,dixs);title('Wave intensity from pressure-only method');
    % xlabel('time(m/s)');ylabel('Wave intensity(W/m^2)')
    %% Comparing velocity with Pxs
    % t = [0:length(velocity)-1]/sampling_rate;
    % t2 = [0:length(Pxs)-1]/sampling_rate;
    % figure; subplot(1,2,1);plot(t,velocity);title('True Velocity waveform');
    % xlabel('time(s)');ylabel('Velocity(m/s)');
    % subplot(1,2,2);plot(t2,Pxs);title('Pxs');
    % xlabel('time(s)');ylabel('Velocity(m/s)');
    
    
    %% Print figures
    % Pressure and central reservoir
%     T_all=1/sampling_rate*(0:(length(P_all)-1));
    % TimeC=(1:length(cP_av))/sampling_rate;
    TimeC=(1:length(P_av))/sampling_rate;
    %clf                  
    f1 = figure(1); set(f1,'Color','w');
    plot(TimeC,P_av,TimeC,Pr_av,'r', TimeC,Pxs,'k');
    xlabel('Time (s)')
    ylabel('BP (mmHg)')
    title('P, Pres, Pxs')
    box off;
    
    
    %% Variable calculations

    [max_p_av,max_tp_av]=max(P_av);
    [max_pxs,max_txs]=max(Pxs); % calculation
    [max_pr,max_tr]=max(Pr_av); % calculation
%     Pr_less_dia=Pr_av-min(Pr_av);
%     [max_prld,max_trld]=max(Pr_less_dia); % calculation
    euc_dist = dtw(di,dixs);        % DTW distance 
    [max_velocity,time_max_velocity] = max(velocity);
    
    
    
    %% write variables
    proc_var{record_no,1}= filename;         % filename
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
    proc_var{record_no,25}=max_velocity;        % Maximum Velocity
    proc_var{record_no,26}=time_max_velocity/sampling_rate;   % Time of maximun velocity
    proc_var{record_no,27}= max_txs/sampling_rate;            % Time of maximun Pxs
    proc_var{record_no,28} = version;            % Version of reservior analysis used

    
    fprintf('Processing............... %f%% \n', file_number/no_of_files.*100)
%     file_number
end


%% Save the results as an excel spreadheet
datafolder='C:\Spdata\results\';
if ~exist(datafolder, 'dir')
    mkdir(datafolder);
end
xlsfile='C:\Spdata\results\resdata.xlsx'; % Or use .xls for compatibility
header = {'re_file' 'dippks1' 'dimpks' 'dippks2' 'dippksxs1' 'dimpksxs'...
    'dippksxs2','dimlocs','dimlocsxs',	're_sam_rate', 're_maxxsp'...
    're_tn' 're_pinf' 're_pn' 'wri' 'wrixs','re_pb_pf', ...
    'diploc1','diploc2','diplocxs1','diplocxs2','corr_Pxs', 'corr_dI', ... 
    'DTW_dist','max_v','t_vmax','txs_max','kres_v'}; % header

% writetable
Results_table=cell2table(proc_var, 'VariableNames',header);
writetable(Results_table, xlsfile);

clc;
fprintf('Analysis Complete')
