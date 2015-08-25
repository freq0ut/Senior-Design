close all;
clear all;
clc;

fS = 100;
tS = 1/fS;

N0 = 2^8;

t = 0:tS:(N0-1)*tS;
y = 10*sin(2*pi*1*t);
y = awgn(y,10);

stem(y)

pks = FIND_LOCAL_MAXES(y,110);