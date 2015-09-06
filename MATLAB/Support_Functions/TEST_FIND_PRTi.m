close all;
clear all;
clc;

fSCE = 30E+3;
PRT = 2;
PW = 1.3E-3;

fADC = 300E+3;
tADC = 1/fADC;

t = 0:tADC:PW;
y = sin(2*pi*fSCE*t);
y = padarray(y,[0,1E+3],0,'pre');
y = padarray(y,[0,round((PRT-PW)/tADC)-1E+3],0,'post');
y = repmat(y,[1,2]);

% THD = 0.01;
% PRTi = FIND_PRTi(y,THD);

figure(1);
    stem(y,'-b');
%     hold on;
%     line([1,length(y)],[THD,THD],'Color',[1,0,0]);
%     hold off;