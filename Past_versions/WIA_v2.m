%% Wave Intensity Analysis with parameters output.
function [di,wri,diplocs,dippks,dimlocs,dimpks,dipt,di_plus,di_minus] = WIA_v2(Npoly,Frame,P,U,sampling_rate,wavespeed)
%%%%    Versions 1.0, working version with peak location/value/area output.
%%%%    Version 2.0, added wave separation using input wavespeed. Also
%%%%    calculate peaks accroding to the separated waves.

% P - Pressure waveform
% U - velocity waveform

% calculate derivatives
[b,g]=sgolay(Npoly,Frame);   % Calculate S-G coefficients
HalfWin=((Frame+1)/2) -1;
% p=P-min(P);
% p = zeros(1,length(P_av));
% p(1:11)=P_av(end-10:end);
% p(12:end)=P_av(1:end-11);
% p=(p*mmHgPa);
% u = zeros(1,length(velocity));
% u(1:11)=U(end-10:end);
% u(12:end)=U(1:end-11);
% u=uconst*(u/max(u));    %Make velocity unitless
% u=u/100;       % retain P and U as the original data
N=length(P);
dp=zeros(N,1); du=dp; ps=dp; us=dp;
for n=(Frame+1)/2:N-(Frame+1)/2
  % Zero-th derivative (smoothing only) 
  ps(n)=dot(g(:,1),P(n-HalfWin:n+HalfWin));
  us(n)=dot(g(:,1),U(n-HalfWin:n+HalfWin));
  % 1st differential
  dp(n)=dot(g(:,2),P(n-HalfWin:n+HalfWin));
  du(n)=dot(g(:,2),U(n-HalfWin:n+HalfWin));
end

du_minus = 1/2.*du - 1/2.*dp/wavespeed;
du_plus = 1/2.*du + 1/2.*dp/wavespeed;
dp_plus = 1/2.*dp + 1/2.*du.*wavespeed;
dp_minus = 1/2.*dp - 1/2.*du.*wavespeed;

di=dp.*du;
% di_plus = dp_plus.*du_plus; 
di_plus = 1/(4.*wavespeed.*1060).*(dp+wavespeed.*du.*1060).^2;
di_plus = di_plus*length(dp)^2/10000;
di_minus = -1/(4.*wavespeed.*1060).*(dp-wavespeed.*du.*1060).^2;
di_minus = di_minus*length(dp)^2/10000;
% di_minus = dp_minus.*du_minus;
di=di*length(dp)^2/10000;             % units fixed - now in W/m2 * 10^4 per cycle^2

minpeak=max(di_plus)/20;             % peaks <5% of Wf1 ignored
% I've left the warning when there are no peaks for now but if it
% needs to be suppressed then 'signal:findpeaks:largeMinPeakHeight' is
% its id
[~,lsys]=min(dp);               % restrict analysis to systole
% lsys=round(lsys)+15;             % round and add 5 samples to give a margin for error for duration of systole
% [dippks(1),diplocs(1), dipw(1)]=findpeaks(di(1:lsys), 'NPeaks',1,'MinPeakHeight',minpeak); % find 1st dI+ peaks (Wf1)
% [dimpks,dimlocs,dimw]=findpeaks(-di(1:lsys), 'NPeaks',1,'MinPeakHeight',0.7*max(-di)); % find one dI- peaks
% [dippks(2),diplocs(2), dipw(2)]=findpeaks(fliplr(di(1:lsys)), 'NPeaks',1,'MinPeakHeight',minpeak); % find 2nd dI+ peaks (Wf2) by flipping the data and running findpeak in the other direction
% diplocs(2)=lsys-diplocs(2)+1; % correct location for the flipping
lsys=round(lsys)+30;             % round and add 5 samples to give a margin for error for duration of systole
[dippks(1),diplocs(1), dipw(1)]=findpeaks(di_plus(1:lsys), 'NPeaks',1,'MinPeakHeight',minpeak); % find 1st dI+ peaks (Wf1)
[temp1,temp2, temp3]=findpeaks(di_plus(1:lsys), 'NPeaks',2,'MinPeakHeight',minpeak); % find 2nd dI+ peaks (Wf1)
if temp1>dippks(1)
    dippks(1) = temp1;diplocs(1)=temp2;dipw(1)=temp3;
end
[dimpks,dimlocs,dimw]=findpeaks(-di_minus(1:lsys), 'NPeaks',1,'MinPeakHeight',0.5*max(-di_minus)); % find one dI- peaks
[dippks(2),diplocs(2), dipw(2)]=findpeaks(flip(di_plus(1:lsys)), 'NPeaks',1,'MinPeakHeight',minpeak); % find 2nd dI+ peaks (Wf2) by flipping the data and running findpeak in the other direction
diplocs(2)=lsys-diplocs(2)+1; % correct location for the flipping
%% check peaks
%      figure; hold on; plot(dixs); plot(diplocs(1),dippks(1),'ko');
%      plot(dimlocs,-dimpks,'ro'); plot(diplocs(2),dippks(2),'ks');
%%
%     dippks=dippks*mmHgPa*length(dp);
%     dimpks=dimpks*mmHgPa*length(dp); // correct error
%     dippks=dippks*length(dp).^2;
%     dimpks=dimpks*length(dp).^2;
dipt=diplocs/sampling_rate;
dimt=dimlocs/sampling_rate;
% 
% % error trap for when Wf2 is unmeasureable
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


end



