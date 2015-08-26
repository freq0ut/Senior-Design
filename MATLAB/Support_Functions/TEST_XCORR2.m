close all;
clear all;
clc;

fS = 100;
tS = 1/fS;

N0 = 2^10;

t=0:tS:(N0-1)*tS;
y1 = sin(2*pi*1*t);
y2 = sin(2*pi*1*(t+0.123));

[XC, XC_Lags] = XCORR2(y1,y2,30);

stem(XC_Lags,XC);
