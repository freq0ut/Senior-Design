close all;
clc;

fS = 10;
tS = 1/fS;

N0 = 512;

t=0:tS:(N0-1)*tS;
y1 = sin(2*pi*t);
y2 = sin(2*pi*(t+pi/8));

%[XC, XC_Lags] = xcorr(y1,y2,'coeff');

figure(1);
    subplot(1,2,1);
        stem(DATA(:,1),DATA(:,2),'-b');
        xlim([-50,50]);
        ylim([-1,1]);
    subplot(1,2,2);
        stem(DATA(1:101,3),DATA(1:101,4),'-r');
        xlim([-50,50]);
        ylim([-1,1]);
        