clear all;
close all;
clc;

MPD = 30;
xNBRS = 10;

tS = 0.01;
N0 = 1024;

t=0:tS:(N0-1)*tS;
y = sin(2*pi*t);

stem(y);
[pks,pkLocs] = PEAKS_FINDER(y,MPD,xNBRS);