clear all;
close all;
clc;

THD   = 8;
MPD   = 3;
xNBRS = 3;

Y = [7  8  9 10  9  8  7  7  8  9 10  9  8  7  6  5  4  5  6  5  5  4];
X = [1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22];

[pks,pkLocs] = FIND_PEAKS(Y,THD,MPD,xNBRS);

stem (X,Y,'-b');
hold on;
stem(pkLocs,pks,'-r');
hold off;