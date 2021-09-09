%% Function to help locate the last cardiac cycle contained in the dataset
%  Where the last cycle is suppose to be the most stable(Reaching periodic state.)

function  [sysloc, dialoc]=sys_dia_finder(signal,sampling_rate) 

% analyse the whole signal
    p=signal';                  % Non-calibrated pressure trace (calibrate later)
%     j=find(p==0);               % gets rid of any zeros which occasionally appear at the end of the trace (probably due to the transfer function-induced shift)   
%     p=p(1:j(1)-1);
    
    
%     p_av=pulse';                % calibrated pressure trace
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
end