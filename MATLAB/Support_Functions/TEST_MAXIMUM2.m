close all;
clear all;
clc;

y = [0 1 2 3 4 5 6 7 8 9 9 9 8 7 6 5 4 3 2 1 10];
x = 1:length(y);

[yMax, xMax] = MAXIMUM2(x,y,x(1),x(end));
