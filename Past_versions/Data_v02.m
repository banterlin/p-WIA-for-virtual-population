clc; clear all;close all;

% load 'left-radial.txt';
% load 'inlet.txt';

[filename,filepath] = uigetfile('*.txt');
folder_name ='C:\Spdata\single\';
% signal = left_radial;
signal = load(fullfile(filepath,filename));

sampling_rate = 1000;

%% Plotting graphs for patientF
% t = signal(:,1);
% plot_title = {'time(s)', 'Area(cm^2)','Velocity(cm/s)', ...
%             'Pressure(Pa/100)','Wave Speed(cm/s)'};
% hf1=figure(20);set(hf1,'Color','w');
% for n = 2:5
%     subplot(4,1,n-1)
%     plot(t,signal(:,n));xlabel('Time');title(plot_title(n));zoom xon;
% end
%% Output Preesure to external TXT File
% Manual Selection of Cradiac Pulse Timing
% Single pulse
% LRBP(:,1) = left_radial(LP,5);
% LRBP(:,2) = inlet(LP,5);
% LRBP(:,3) = left_radial(LP,4);
% LRBP(:,4) = inlet(LP,4);
% dlmwrite('LRBP.txt',LRBP,'delimiter',' ')

% Finding the last cycle
[sysloc, dialoc]=sys_dia_finder(signal(:,4),sampling_rate) ;
periph_T  = [dialoc(end-1) : dialoc(end)];



% [c_sysloc, c_dialoc]=sys_dia_finder(inlet(:,4),sampling_rate) ;
% central_T  = [c_dialoc(length(c_dialoc)-1) : c_dialoc(length(c_dialoc))];

% Whole Pulse 
LRBP = NaN(length(signal),3);
LRBP(periph_T,1) = signal(periph_T,3); % Velocity profile for central pressure
% LRBP(periph_T,1) = flip(signal(periph_T,3)); %Fliped Velocity profile for central pressure
LRBP(:,2) = signal(:,4);             % Central bp signal(5 cardiac cycles)
LRBP(periph_T,3) = signal(periph_T,4);      % Periphral BP (single cycle)
% LRBP(periph_T,3) = flip(signal(periph_T,4));      % FlipedPeriphral BP (single cycle)



hf2=figure(122);set(hf2,'Color','w');plot(LRBP(periph_T,3));xlabel('Time');title('Pulse of pressure');
% hf3=figure(133);set(hf3,'Color','w');plot(LRBP(periph_T,1));xlabel('Time');title('Pulse of velocity');

% dlmwrite([folder_name filename],LRBP,'delimiter',' ');
dlmwrite([folder_name 'LBRP.txt'],LRBP,'delimiter',' ');
fprintf('Data Extraction done!');


%% Test Saved data
%     fid = fopen('LRBP.txt');
%     data = textscan(fid,'%f%f%f%f%f');
% % Plotting 
% plot(central_T/1000,LRBP(central_T,4))
% 
%  plot(periph_T/1000,LRBP(periph_T,3))
% plot(t, LRBP(:,2));
% fclose(fid);




