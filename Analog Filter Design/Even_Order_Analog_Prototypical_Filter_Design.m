% Author: Joshua Simmons
%
% Date: October 22, 2015
%
% Desciption: Plots the frequency response of even N-th order
%               bandpass filters.

clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% USER INPUT %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fMid  = 20E+3;% Mid frequency
Q = 10;% Quality factor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% CALCULATING TRANSFER FUNCTION %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = linspace(1,2*fMid,1E+4);
s = 1i*2*pi*f;

omega = 2*pi*fMid;

zeta = 1/(2*Q);

% Prototypical 2nd Order LPF Transfer Function
H_LPF2 = (1/Q) ./ (s.^2+2*zeta*omega.*s+omega^2);

% Prototypical 2nd Order HPF Transfer Function
H_HPF2 = (1/Q) * (s.^2) ./ (s.^2+2*zeta*omega.*s+omega^2);

% Prototypical 2nd Order BPF Transfer Function
H_BPF2 = (1/Q) * (omega^2) ./ (s.^2+2*zeta*omega.*s+omega^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% PLOTTING TRANSFER FUNCTION %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plotting Low Pass Filters
figure(1);
    subplot(2,1,1);
        plot(0.001*f,abs(H_LPF2));% 2nd Order
        hold on;
        grid on;
        plot(0.001*f,abs(H_LPF2.^2));%  4th Order
        plot(0.001*f,abs(H_LPF2.^3));%  6th Order
        plot(0.001*f,abs(H_LPF2.^4));%  8th Order
        plot(0.001*f,abs(H_LPF2.^5));% 10th Order
        %ylim([0,1]);
        xlabel('Frequency [kHz]');
        ylabel('|| H(j\omega) ||');
        legend({'2nd','4th','6th','8th','10th'});
        titleString = sprintf('Lowpass Filter Frequency Response\nQ = %3.1f', Q);
        title(titleString);
        hold off;
	subplot(2,1,2);
        plot(0.001*f,(180/pi)*angle(H_LPF2));% 2nd Order
        hold on;
        plot(0.001*f,(180/pi)*angle(H_LPF2.^2));%  4th Order
        plot(0.001*f,(180/pi)*angle(H_LPF2.^3));%  6th Order
        plot(0.001*f,(180/pi)*angle(H_LPF2.^4));%  8th Order
        plot(0.001*f,(180/pi)*angle(H_LPF2.^5));% 10th Order
        ylim([-180,180]);
        set(gca,'YTick',[-180,-135,-90,-45,0,45,90,135,180])
        set(gca,'YTickLabel',[-180,-135,-90,-45,0,45,90,135,180])
        grid on;
        xlabel('Frequency [kHz]');
        ylabel('\Theta H(j\omega) [Deg]');
        legend({'2nd','4th','6th','8th','10th'});
        hold off;

% Plotting High Pass Filters
figure(2);
    subplot(2,1,1);
        plot(0.001*f,abs(H_HPF2));% 2nd Order
        hold on;
        grid on;
        plot(0.001*f,abs(H_HPF2.^2));%  4th Order
        plot(0.001*f,abs(H_HPF2.^3));%  6th Order
        plot(0.001*f,abs(H_HPF2.^4));%  8th Order
        plot(0.001*f,abs(H_HPF2.^5));% 10th Order
        ylim([0,1]);
        xlabel('Frequency [kHz]');
        ylabel('|| H(j\omega) ||');
        legend({'2nd','4th','6th','8th','10th'});
        titleString = sprintf('Highpass Filter Frequency Response\nQ = %3.1f', Q);
        title(titleString);
        hold off;
	subplot(2,1,2);
        plot(0.001*f,(180/pi)*angle(H_HPF2));% 2nd Order
        hold on;
        plot(0.001*f,(180/pi)*angle(H_HPF2.^2));%  4th Order
        plot(0.001*f,(180/pi)*angle(H_HPF2.^3));%  6th Order
        plot(0.001*f,(180/pi)*angle(H_HPF2.^4));%  8th Order
        plot(0.001*f,(180/pi)*angle(H_HPF2.^5));% 10th Order
        ylim([-180,180]);
        set(gca,'YTick',[-180,-135,-90,-45,0,45,90,135,180])
        set(gca,'YTickLabel',[-180,-135,-90,-45,0,45,90,135,180])
        grid on;
        xlabel('Frequency [kHz]');
        ylabel('\Theta H(j\omega) [Deg]');
        legend({'2nd','4th','6th','8th','10th'});
        hold off;

% Plotting Band Pass Filters
figure(3);
    subplot(2,1,1);
        plot(0.001*f,abs(H_BPF2));% 2nd Order
        hold on;
        grid on;
        plot(0.001*f,abs(H_BPF2.^2));%  4th Order
        plot(0.001*f,abs(H_BPF2.^3));%  6th Order
        plot(0.001*f,abs(H_BPF2.^4));%  8th Order
        plot(0.001*f,abs(H_BPF2.^5));% 10th Order
        ylim([0,1]);
        xlabel('Frequency [kHz]');
        ylabel('|| H(j\omega) ||');
        legend({'2nd','4th','6th','8th','10th'});
        titleString = sprintf('Bandpass Filter Frequency Response\nQ = %3.1f', Q);
        title(titleString);
        hold off;
	subplot(2,1,2);
        plot(0.001*f,(180/pi)*angle(H_BPF2));% 2nd Order
        hold on;
        plot(0.001*f,(180/pi)*angle(H_BPF2.^2));%  4th Order
        plot(0.001*f,(180/pi)*angle(H_BPF2.^3));%  6th Order
        plot(0.001*f,(180/pi)*angle(H_BPF2.^4));%  8th Order
        plot(0.001*f,(180/pi)*angle(H_BPF2.^5));% 10th Order
        ylim([-180,180]);
        set(gca,'YTick',[-180,-135,-90,-45,0,45,90,135,180])
        set(gca,'YTickLabel',[-180,-135,-90,-45,0,45,90,135,180])
        grid on;
        xlabel('Frequency [kHz]');
        ylabel('\Theta H(j\omega) [Deg]');
        legend({'2nd','4th','6th','8th','10th'});
        hold off;
        