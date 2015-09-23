close all;
clear all;
clc;

X = 0:15;
Y = 0:15;

A = SIMPSON(Y,X)

stem(X,Y)