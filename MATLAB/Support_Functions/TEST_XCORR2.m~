close all;
clear all;
clc;

fS = 100;
tS = 1/fS;

N0 = 2^10;

t=0:tS:(N0-1)*tS;
y = sin(2*pi*1*t);

iBreak = BREAK_THRESHOLD(y,0.8);
tBreak = iBreak*tS;

plot(y);