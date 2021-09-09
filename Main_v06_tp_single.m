%% Script adapted from main_v06 to fit true patients instead of virtual ones
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
%  v0.6 Added wave separation using true wavespeed/estimated wavespeed from
%  1ms-1 assumption. Added Diastolic Excess pressure zeroing function.
%  Added number of zero crossing on the dP as parameter.

%%%%%%%%%%%%%%%
%% m files required to be in directory
% fitres_v6.m
% kreservoir_v14.m
%%%%%%%%%%%%%%%%
%% Constants
clear all; clc;

close all;


sampling_rate = 128;  %

kres_v=14;          % Version tracking for reservoir fitting
headernumber=25;       % headers for columns of results (see end)
mmHgPa = 133;          % P conversion for WIA
uconst=1;              % empirical constant to convert normalized
% velocity to m/s based on Hughes et al.
% Front. Physiol. DOI: 10.3389/fphys.2020.00550
Npoly=3;               % Order of polynomial fit for sgolay
Frame=11;               % Window length for sgolay
version='0.6';         % Version of Main
%%%%%%%%%%%%%%%%
%% Select files
% default folder as per manual
folder_name ='C:\Spdata\Real\single\';
% check that C:\Spdata exists and if not allows new folder to be chosen
if ~exist('C:\Spdata\Real\single\', 'dir')
    answer = questdlg('C:\Spdata\Real doesnt exist. Would you like to choose another folder?', ...
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
for file_number=1:record_no
    % refresh filename
%     record_no = file_number;        %for multiple file input
    filename=file_lists(record_no).name;
    
    fid = fopen([folder_name filename]);
%     data = textscan(fid,'%f%f%f%f%f','headerlines',5);
    data = dlmread([folder_name filename],'\t',4,0);
    fclose(fid);

    pressure=data(:,3);
%     index = find(pressure == 0, 1, 'first');
%     pressure = pressure(1 : index-1);
    
    %     central_pulse=data{1,4};
    velocity = data(:,4)./100 ; % Convert cm/s in raw data to m/s
% Following code checks the peaks between 1 detected cycle with 1.5
% detected cycle, to avoid cases where young subjects showed very low
% magnitude of dirotic notch (where the algorithm would see it as a diastolic trough)
%     [sysloc, dialoc]=sys_dia_finder(pressure,sampling_rate) ;
%     if max(pressure(dialoc(end-2) : dialoc(end-1))) - max(pressure(dialoc(end-3) : dialoc(end-2))) > 20
%         periph_T  = [dialoc(end-4) : dialoc(end-2)];
%         fprintf('1');
%     elseif max(pressure(dialoc(end-2) : dialoc(end-1))) - max(pressure(dialoc(end-3) : dialoc(end-2))) < -20
%         periph_T  = [dialoc(end-3) : dialoc(end-1)];
%         fprintf('2');
%     else
%         periph_T  = [dialoc(end-2) : dialoc(end-1)];
%         fprintf('3');
%     end
    
%      periph_T  = [dialoc(end-1) : dialoc(end)];

    
%     single_pulse = pressure(periph_T,1);    % Single cardiac cycle analysis
%     velocity = velocity(periph_T,1);        % Single cardiac cycle analysis
    
%     single_pulse=single_pulse(~isnan(single_pulse));     % remove NaNs
    single_pulse=data(:,3);
%     central_signal= central_signal(~isnan(central_signal));
    velocity = velocity(1:length(single_pulse));
    %     central_pulse=central_pulse(~isnan(central_pulse));     % remove NaNs
    
    % call function for reservoir fitting
    [P_av, Pr_av,Pinf_av,Pn_av,Tn_av,Sn_av, fita_av, fitb_av,rsq_av,dp,nn]=fitres_v6(single_pulse,sampling_rate,kres_v);
    
    Pxs=P_av-Pr_av; % Excess Pressure
    % If Pxs < Threshold , we use delayed fitting time (V16)
%     if max(Pxs(1:400))- max(Pxs(400:end)) < 5; % Threshold Value here require more study
%         kres_v = 16;
%         [P_av, Pr_av,Pinf_av,Pn_av,Tn_av,Sn_av, fita_av, fitb_av,rsq_av,dp,nn]=fitres_v6(single_pulse,sampling_rate,kres_v);
%         Pxs=P_av-Pr_av; % Excess Pressure
%     end
    zci = @(v) find(diff(sign(v))); % Function for finding zero crossing
    n_Zx = length(zci(dp));         % Number of zero crossing presented by the pressure wave
    
%     f99 = figure(99); set(f99,'Color','w');plot(dp);title('dp')% Figure for rate of change of pressure
    
    %% Diastolic Excess pressure zeroing. 
%     idx1= zci(Pxs);
%     if length(idx1) > 2 ;
%         idx2 = idx1(3);
%     else
%         idx2 = idx1(2);
%     end
%     Pxs(idx2:end) = 0;
    %% PU loop
%     f88 = figure(88); set(f88,'Color','w');
%     plot(Pxs,P_av) ; title('PULoop for Pxs')
%     grid on;
%     f89 = figure(89); set(f89,'Color','w');
%     plot(velocity(1:length(P_av)),P_av*mmHgPa/1000) ; title('PULoop')
%     grid on;xlabel('ms^{-1}');ylabel('kPa');
    
 
    
    %% Print Figure for reservior analysis
    
    TimeC=(1:length(P_av))/sampling_rate;

    % figure('visible','off');                     % dont display figure
    f1 = figure(1); set(f1,'Color','w');

    plot(TimeC,P_av,TimeC,Pr_av,'r', TimeC,Pxs,'k');legend('P','Pres','Pxs')
    xlabel('Time (s)')
    ylabel('BP (mmHg)')
    title('P, Pres, Pxs')
    box off;
    
    % calculate mean arterial pressure from waveform
    %map=str2double(spdata{2,32}); % data in file
    cPb_av=(Pr_av-min(P_av))/2;
    cPf_av=P_av-cPb_av-min(P_av);
    Pb_Pf=max(cPb_av)/max(cPf_av);
    RI=max(cPb_av)/(max(cPf_av)+ max(cPb_av));
    
    % duration of diastole
    % Tdia=(length(P_av)-Sn_av)/sampling_rate;
    
    
    
    
    %% Wave intensity analysis (WI Using pressure-only method)
    cpwia = zeros(1,length(P_av));
%     cpwia(1:11)=P_av(end-10:end);
%     cpwia(12:end)=P_av(1:end-11);

    cpwia = P_av;
    
    cpwia=(cpwia*mmHgPa);
    cuxwia = zeros(1,length(Pxs));
%     cuxwia(1:11)=Pxs(end-10:end);
%     cuxwia(12:end)=Pxs(1:end-11);
    cuxwia = Pxs;    
    % convert xwia to flow velocity and assume peak velocity = 1m/s
    % based on Lindroos, M., et al. J Am Coll Cardiol 1993; 21: 1220-1225.
    % cuxwia=uconst*(cuxwia-min(cuxwia))/(max(cuxwia)-min(cuxwia));
    
%     uconst = max(velocity); % Using actual max(U) for better estimation
    cuxwia=uconst*(cuxwia/max(cuxwia));   % Peak Velocity calibrated to 1m/s assumption/max(U) depend on uconst
    % Estimate c (wavespeed) as k*dP/du where k is empirical constant
    % currently k = 1!
    rhocxs=max(Pxs)*mmHgPa./1060./uconst; % fixed units (kg) 
    [dixs,wrixs,diplocsxs,dippksxs,dimlocsxs,dimpksxs,diptxs,dimtxs,dipareaxs,dimareaxs] = WIA(Npoly,Frame,cpwia,cuxwia,sampling_rate);
   
    
    
    %% Corrcoef between Pxs and velocity from raw data.

    fprintf("r^2 of diastolic fit is %f \n",rsq_av);
    % Comparing wave intensities
    % Constructing time interval over sampling period
    t = [0:length(dixs)-1]./sampling_rate;
    f2 = figure(2);set(f2,'Color','w'); 
    plot(TimeC,dixs);title('Wave intensity from pressure-only method');
    xlabel('time(m/s)');ylabel('Wave intensity(W/m^2)')
    
    
    
    %% Variable calculations

    [max_p_av,max_tp_av]=max(P_av);
    [max_pxs,max_txs]=max(Pxs); % calculation
    [max_pr,max_tr]=max(Pr_av); % calculation
%     Pr_less_dia=Pr_av-min(Pr_av);
%     [max_prld,max_trld]=max(Pr_less_dia); % calculation
    max_velocity = max(velocity);
    

end