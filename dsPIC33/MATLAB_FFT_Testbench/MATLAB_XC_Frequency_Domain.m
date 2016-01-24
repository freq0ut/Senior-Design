close all;
clear all;
clc;

N = 1024;

fS = 1e6;
Ts = 1/fS;

f0 = fS/(N-1);

t = Ts*(0:N-1);
f = f0*(-N/2:N/2-1);

y1 = cos(2*pi*40e3*t);
y2 = cos(2*pi*40e3*(t+5e-6));

y2R = fliplr(y2);

XC_f = fft(y1,2*N) .* fft(y2R,2*N);

XC_t = ifft(XC_f);
XC_t = XC_t ./ max(abs(XC_t));

XC_t(end) = [];

fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(3,1,1);
        plot(Ts*(-N+1:N-1),XC_t,'-b');
    subplot(3,1,2);
        plot(Ts*(-N+1:N-1),xcorr(y1,y2,'coef'),'-r');
    subplot(3,1,3);
        plot(Ts*(-N+1:N-1),XC_t-xcorr(y1,y2,'coef'),'-g');