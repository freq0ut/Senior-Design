% Joshua Simmons
%
% January 11, 2016
%
% Crosschecking Microchip DSPic FFT Output with a square wave input.

close all;
clear all;
clc;

N = 256;

fS = 1e3;
Ts = 1/fS;

f0 = 1/((N-1)*Ts);

f = f0*(-N/2:N/2-1);
t = Ts*(0:N-1);

y = ones(1,N);
for i=0:N-1;
    if (mod(i,10) >= 5)
        y(i+1) = -1;
    end
end

Y = fftshift(fft(y)) / N;

fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(2,2,1);
    	stem(t,y,'-b');
    	grid on;
    	grid minor;
        xlim([t(1),t(end)]);
    	xlabel('Time [s]');
    	ylabel('Amplitude [V]');
        titleString = sprintf('TEMPORAL\nMATLAB Simulated Input Signal');
        title(titleString);
    subplot(2,2,3);
    	stem(f,abs(Y),'-r');
    	grid on;
    	grid minor;
    	xlabel('Frequency [Hz]');
    	ylabel('Magnitude [Volts]');
        titleString = sprintf('FREQUENCY\nMATLAB Simulated Input Signal');
        title(titleString);