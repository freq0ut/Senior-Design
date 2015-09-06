close all;
clear all;
clc;

N0 = 512;

fADC = 1800E+3;
tADC = 1/fADC;

fSCE = 30E+3;
tSCE = 1/fSCE;

t = 0:tADC:(N0-1)*tADC;
y = sin(2*pi*fSCE*t);
y = padarray(y,[0,1E+2],0,'pre');
y = padarray(y,[0,1E+3],0,'post');

THD = 0.5;
iBreak = BREAK_THRESHOLD(y,THD,'RL');

figure(1);
    stem(y);
    hold on;
    line([1,length(y)],[THD,THD],'Color',[1,0,0]);
    hold off;