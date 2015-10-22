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

fMid  = 40E+3;% Mid frequency
Q = 10;% Quality factor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% CALCULATING TRANSFER FUNCTION %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = linspace(1,2*fMid,1E+3);
s = 1i*2*pi*f;

omega = 2*pi*fMid;

A = (1/Q) * (omega^2) ./ (s.^2+(omega/Q).*s+omega^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% PLOTTING TRANSFER FUNCTION %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1);
    subplot(2,1,1);
        plot(0.001*f,abs(A));% 2nd Order
        hold on;
        grid on;
        plot(0.001*f,abs(A.^2));%  4th Order
        plot(0.001*f,abs(A.^3));%  6th Order
        plot(0.001*f,abs(A.^4));%  8th Order
        plot(0.001*f,abs(A.^5));% 10th Order
        ylim([0,1]);
        xlabel('Frequency [kHz]');
        ylabel('|| H(j\omega) ||');
        legend({'2nd','4th','6th','8th','10th'});
        titleString = sprintf('Bandpass Filter Frequency Response\nQ = %3.1f', Q);
        title(titleString);
        hold off;
	subplot(2,1,2);
        plot(0.001*f,(180/pi)*angle(A));% 2nd Order
        hold on;
        plot(0.001*f,(180/pi)*angle(A.^2));%  4th Order
        plot(0.001*f,(180/pi)*angle(A.^3));%  6th Order
        plot(0.001*f,(180/pi)*angle(A.^4));%  8th Order
        plot(0.001*f,(180/pi)*angle(A.^5));% 10th Order
        ylim([-180,180]);
        grid on;
        xlabel('Frequency [kHz]');
        ylabel('\Theta H(j\omega) [Deg]');
        legend({'2nd','4th','6th','8th','10th'});
        hold off;
        