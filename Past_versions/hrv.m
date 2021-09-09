%% hrv_v1 - batch analysis of radial pulse data using kreservoir_vXX by KHP
%  Copyright 2019 Alun Hughes
%% Modified by Lin
%  This software is distributed under under the terms of the GNU General Public License
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
% http://www.gnu.org/licenses/gpl.html

%% Versions
% v1 First stable version(11/01/20)
%%
function  [P_all,sdsbp_mmhg ,sysloc, dialoc]=hrv(pulse,signal,sampling_rate) 

% analyse the whole signal
    p=signal';                  % Non-calibrated pressure trace (calibrate later)
%     j=find(p==0);               % gets rid of any zeros which occasionally appear at the end of the trace (probably due to the transfer function-induced shift)   
%     p=p(1:j(1)-1);
    
    
    p_av=pulse';                % calibrated pressure trace
    % L=length(p);
    pnorm=(p-min(p))/(max(p)-min(p));
    t=(0:length(p)-1)/sampling_rate;    
    
%     % find systolic peaks
%     [syspeaks, sysloc]=findpeaks(pnorm,'MinPeakDistance',64); % enforces that HR < 120bpm
%     % find diastolic peaks
%     [~, dialoc]=findpeaks(abs(pnorm-1),'MinPeakDistance',64);
%     diapeaks=pnorm(dialoc);
    promfact=0.08; % was 0.25
    [syspeaks, sysloc]=findpeaks(p, 'MinPeakProminence',max(p)*promfact,'MinPeakDistance',50);
    promfact=0.25;
    p_upside=-p-min(-p);
    [~, dialoc] = findpeaks(p_upside, 'MinPeakProminence',max(p_upside)*promfact,'MinPeakDistance',30);
    clear p_upside;
    diapeaks=pnorm(dialoc);
   
    % calibrate p
    avsbp=max(p_av); avdbp=min(p_av); avpp= avsbp-avdbp;
    avSsig=mean(syspeaks); avDsig=mean(diapeaks); avPPsig=avSsig-avDsig;
    P_all=((pnorm-mean(diapeaks))*avpp)+avdbp;
    sbp=P_all(sysloc);
    dbp=P_all(dialoc);
    sdsbp_mmhg= std(sbp(1:end-1)); % Standard deviation of SBP, mmHg   

end