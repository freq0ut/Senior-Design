clear all;
close all;
clc;

N = 2^10;

fS = 1E+6;
tS = 1/fS;

t = 0:tS:(N-1)*tS;

chan1 = sin(2*pi*40000*t);
chan2 = sin(2*pi*40000*(t-20E-6));

%iStart = -N+1;
%iStop = N-1;

iStart = floor(-30E-6/tS);
iStop = floor(-10E-6/tS);

tic
[myXC, myXC_Lags] = XCORR3(chan1,chan2,iStart,iStop);
toc

tic
[MATLAB_XC, XC_Lags] = xcorr(chan1,chan2);
toc

error = MATLAB_XC - myXC;

figure(1);
    subplot(2,2,1);
        stem(myXC_Lags,myXC,'-b');
    subplot(2,2,3);
        stem(XC_Lags,MATLAB_XC,'-r');
    subplot(1,2,2);
        plot(XC_Lags,error,'-g');