close all;
clear all;
clc;

addpath('/Users/betio32/Desktop/Senior Design/Senior-Design/MATLAB/Support_Functions');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% INITIALIZATION OF PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Global Simulation Parameters
trialTotal = 10; % Total number of iterations of main loop
DATA_AZ = zeros(2,trialTotal);

% Pinger Properties
SNR  = 20;         % Signal to Noise Ratio [dB]
fPing = 30e3;      % Source freq [Hz]
tPing = 1/fPing;   % Source period [s]
vP   = 1482;       % Propagation Velocity [m/s]
lambda = vP/fPing; % Wavelength [m]
pingMaxDist = 1;   % Pinger max distance from sensors [m]

% Hydrophone Properties
D = lambda/3;    % Hydrophone spacing [m]
d = D/sqrt(2);   % For coordinates of hydrophones [m]

% ADC
fADC = 1.8e6;  % Sample freq [Hz]
tADC = 1/fADC; % Sample period [s]
N0 = 2^11;     % Samples per frame

% Microcontroller Properties
tD_Act = [0;0;0;0]; % Actual time delays
tD_Est = [0;0;0;0]; % Estimated time delays (Trapezoidal Rule)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% CONTRUCTING INPUT SIGNALS %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% START MAIN LOOP
for trialCount = 1:trialTotal;
    % Pinger location (xP,yP,zP)
    Ping_Act(1) = pingMaxDist*(2*rand()-1);
    Ping_Act(2) = pingMaxDist*(2*rand()-1);
    Ping_Act(3) = pingMaxDist*(2*rand()-1);

    % Calculating actual azimuths to source
    azimuthH_Act = wrapTo2Pi(atan2(Ping_Act(2),Ping_Act(1))) * (180/pi);
    azimuthV_Act = wrapTo2Pi(atan2(Ping_Act(3),Ping_Act(1))) * (180/pi);
    
    % Calculating the actual sphere radii
    R_Act(1) = sqrt( (Ping_Act(1)-d)^2 + (Ping_Act(2)  )^2 + (Ping_Act(3)  )^2 );
    R_Act(2) = sqrt( (Ping_Act(1)  )^2 + (Ping_Act(2)-d)^2 + (Ping_Act(3)  )^2 );
    R_Act(3) = sqrt( (Ping_Act(1)  )^2 + (Ping_Act(2)  )^2 + (Ping_Act(3)-d)^2 );
    R_Act(4) = sqrt( (Ping_Act(1)-d)^2 + (Ping_Act(2)-d)^2 + (Ping_Act(3)-d)^2 );

    % Determining the actual time delays
    TOA_Act = R_Act(1)/vP;
    tD_Act(2) = (R_Act(2)-R_Act(1)) / vP;
    tD_Act(3) = (R_Act(3)-R_Act(1)) / vP;
    tD_Act(4) = (R_Act(4)-R_Act(1)) / vP;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN SIGNAL PROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Constructing estimated time delays
    error = tADC;
    
    for i=2:4;
        tD_Est(i) = tD_Act(i) + (2*rand()-1) * error;
    end

    % Calculating the estimated Time-Of-Arrival
    TOA_Est  = FIND_TETRA_TOA(d,tD_Est(2),tD_Est(3),tD_Est(4),vP);
    
    % Calculating estimated sphere radii
    R_Est(1) = vP*TOA_Est;
    R_Est(2) = vP*(TOA_Est+tD_Est(2));
    R_Est(3) = vP*(TOA_Est+tD_Est(3));
    R_Est(4) = vP*(TOA_Est+tD_Est(4));
    
    % Determining the estimated source location
    Ping_Est(1) = (-R_Est(1)^2 + R_Est(2)^2 + R_Est(3)^2 - R_Est(4)^2 + 2*d^2) / (4*d);
    Ping_Est(2) = ( R_Est(1)^2 - R_Est(2)^2 + R_Est(3)^2 - R_Est(4)^2 + 2*d^2) / (4*d);
    Ping_Est(3) = ( R_Est(1)^2 + R_Est(2)^2 - R_Est(3)^2 - R_Est(4)^2 + 2*d^2) / (4*d);
    
    % Calculating estimated azimuths to source
    azimuthH_Est = wrapTo2Pi(atan2(Ping_Est(2),Ping_Est(1))) * (180/pi);
    azimuthV_Est = wrapTo2Pi(atan2(Ping_Est(3),Ping_Est(1))) * (180/pi);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROCESSING DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    DATA_AZ(1,trialCount) = azimuthH_Act - azimuthH_Est;
    DATA_AZ(2,trialCount) = azimuthV_Act - azimuthV_Est;
end

string1 = sprintf('f_{ADC} = %0.1f [MHz]', fADC/1E+6);
string2 = sprintf('D = \\lambda/3');

figure(1)
    subplot(1,2,1);
        boxplot(DATA_AZ(1,:),'labels',{'Horizontal'});
        hold on;
        line([-1,1],[5,5]);
        line([-1,1],[-5,-5]);
        ylabel('Azimuth Error (deg)');
        legend({string1,string2});
        hold off;
    subplot(1,2,2);
        boxplot(DATA_AZ(2,:),'labels',{'Vertical'});
        hold on;
        line([-1,1],[5,5]);
        line([-1,1],[-5,-5]);
        legend({string1,string2});
        hold off;

% 5 degrees
% D = 0.10,     fADC =  25.0e6 mediocre
% D = 0.40,     fADC =   1.8e6 more mediocre
% D = lambda/3, fADC = 800.0e6