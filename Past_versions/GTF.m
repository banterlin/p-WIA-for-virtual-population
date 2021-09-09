%% Transfer Function Using ARX model to estimate cBP from pBP %%%%%%%%

%% Version 1.0
clear all; clc; close all;

%constant
mmHgPa = 133;
sampling_rate = 1000;

%% Obtain single pulse of cBP and pBP
folder_name ='C:\Users\xq633\Desktop\Msc\Research Project\P-only_wia\virtualPopulation2\M\40-49\4\';

filename = 'inlet.txt';
fid = fopen([folder_name filename]);
data = textscan(fid,'%f%f%f%f%f','headerlines',5);
fclose(fid);

cBP = data{1,4}/mmHgPa*100;
    [sysloc, dialoc]=sys_dia_finder(cBP,sampling_rate) ;
    if max(cBP(dialoc(end-2) : dialoc(end-1))) - max(cBP(dialoc(end-3) : dialoc(end-2))) >5
        periph_T  = [dialoc(end-4) : dialoc(end-2)];
    elseif max(cBP(dialoc(end-2) : dialoc(end-1))) - max(cBP(dialoc(end-3) : dialoc(end-2))) < -5
        periph_T  = [dialoc(end-3) : dialoc(end-1)];
    else
        periph_T  = [dialoc(end-2) : dialoc(end-1)];
    end
cBP = cBP(periph_T);

filename = 'left-radial.txt';
fid = fopen([folder_name filename]);
data = textscan(fid,'%f%f%f%f%f','headerlines',5);
fclose(fid);

pBP = data{1,4}/mmHgPa*100;

    [sysloc, dialoc]=sys_dia_finder(pBP,sampling_rate) ;
    if max(pBP(dialoc(end-2) : dialoc(end-1))) - max(pBP(dialoc(end-3) : dialoc(end-2))) >5
        periph_T  = [dialoc(end-4) : dialoc(end-2)];
    elseif max(pBP(dialoc(end-2) : dialoc(end-1))) - max(pBP(dialoc(end-3) : dialoc(end-2))) < -5
        periph_T  = [dialoc(end-3) : dialoc(end-1)];
    else
        periph_T  = [dialoc(end-2) : dialoc(end-1)];
    end

pBP = pBP(periph_T);


%% transfer function modelling

A = [1  -1.5  0.7];
B = [0 1 0.5];
sys0 = idpoly(A,B);

y = sim(sys0,[cBP pBP]);

z = [y,cBP];
sys = arx(z,[2 2 1]);


 mcBP = tf(sys);   % transform to tf



