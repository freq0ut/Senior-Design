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
N0 = 2^10;      % Samples per frame

% Microcontroller Properties
MPD = 60;  % Minimum Peak Distance
xNBRS = 2; % Neighbors to look to the left and right of
THD = 0.5; % Threshold

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTRUCTING TEST SIGNAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = 0:tADC:(N0-1)*tADC;
y = cos(2*pi*fPing*t);

iStart = round(N0*rand() - fADC/PRF);

if iStart <= 0
    iStart = 1;
end

for i = iStart:iStart+round(fADC/PRF);
    y(i) = 0;
end

y = awgn(y,SNR);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ESTIMATING PRF PERIOD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[pks, pkLocs] = FIND_PEAKS(y,MPD,xNBRS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VISUALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1);
    stem(t,abs(y));
    grid on;
    xlabel('Time [us]');
    ylabel('Aplitude [V]');

