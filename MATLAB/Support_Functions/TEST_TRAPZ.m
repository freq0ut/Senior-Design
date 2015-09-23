close all;
clear all;
clc;

N0 = 1024;

fADC = 100;

t = 0:1/fADC:(N0-1)/fADC;
y = sin(2*pi*t);

stem(t,y)