%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:  Joshua Simmons
% Started: August, 2015
% Status:  Incomplete
%
% Description: This MATLAB script calculates the PRF Period of the pinger.
%
%              The ADC will sample at a lower frequency to consume less
%              memory. Once the PRF Period is determined, the ADC can
%              step up its sampling frequency and start sampling just
%              before the pinger transistions to ON.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;
clc;

addpath('C:\Users\Joshua Simmons\Desktop\Senior_Design\Senior-Design\MATLAB\Support_Functions');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETER INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Global Simulation Parameters

% Pinger Properties
SNR  = 20;         % Signal to Noise Ratio [dB]
fPing = 30E+3;     % Source freq [Hz]
tPing = 1/fPing;   % Source period [s]
PRF = 1E+3;        % Pulse-Repetitive-Frequency [Hz]
PRT = 1/PRF;       % Pulse-Repetitive-Period [s]
vP = 1482;         % Propagation Velocity [m/s]
lambda = vP/fPing; % Wavelength [m]

% ADC
fADC = 300E+3; % Sample freq [Hz]
tADC = 1/fADC;  % Sample period [s]
N0 = 2^9;      % Samples per frame

% Microcontroller Properties
MPD = 1;  % Minimum Peak Distance
xNBRS = 1; % Neighbors to look to the left and right of
THD = 0.5; % Threshold

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTRUCTING TEST SIGNAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = tADC:tADC:(N0)*tADC;    % Making time start at zero threw off pkLocs by 1
y = cos(2*pi*fPing*t);

iGap = round(fADC/PRF);
iStart = round((N0-iGap)*rand());

for i = iStart:iStart+iGap;
    y(i) = 0;
end

y = awgn(y,SNR);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ESTIMATING PRF PERIOD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[pks, pkLocs] = FIND_PEAKS(abs(y),THD,MPD,xNBRS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VISUALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1);
    subplot(1,2,1);
        stem(t,abs(y));
        hold on;
        stem(pkLocs*tADC,pks,'-r');
        line([t(1), t(end)],[THD, THD],'Color',[0.5,0,0.5]);
        grid on;
        xlabel('Time [us]');
        ylabel('Aplitude [V]');
        hold off;
    subplot(1,2,2);
        stem(abs(y));
        hold on;
        stem(pkLocs,pks,'-r');
        line([1, length(y)],[THD, THD],'Color',[0.5,0,0.5]);
        grid on;
        xlabel('Time [us]');
        ylabel('Aplitude [V]');
        hold off;