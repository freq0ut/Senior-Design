close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% INITIALIZATION OF PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Global Simulation Parameters
trialTotal = 1e4; % Total number of iterations of main loop
DATA_AZ = zeros(1,trialTotal);

% Pinger Properties
fPing = 47e3;      % Source freq [Hz]
tPing = 1/fPing;   % Source period [s]
vP   = 1482;       % Propagation Velocity [m/s]
lambda = vP/fPing; % Wavelength [m]
pingMaxDist = 1;   % Pinger max distance from sensors [m]

% Hydrophone Properties
D = 0.05;    % Hydrophone spacing [m]

% ADC
fADC = 300E+3;  % Sample freq [Hz]
tADC = 1/fADC;  % Sample period [s]
N0 = 1024;      % Samples per frame

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
    
    % Determining the actual time delays    
    R_Act(1) = sqrt( (Ping_Act(1)  )^2 + (Ping_Act(2)  )^2 + (Ping_Act(3)  )^2 );
    R_Act(2) = sqrt( (Ping_Act(1)  )^2 + (Ping_Act(2)+D)^2 + (Ping_Act(3)  )^2 );
    R_Act(3) = sqrt( (Ping_Act(1)+D)^2 + (Ping_Act(2)+D)^2 + (Ping_Act(3)  )^2 );
    R_Act(4) = sqrt( (Ping_Act(1)+D)^2 + (Ping_Act(2)  )^2 + (Ping_Act(3)  )^2 );
    
    TOA_Act = R_Act(1)/vP;
    tD_Act(2) = (R_Act(2)-R_Act(1)) / vP;
    tD_Act(3) = (R_Act(3)-R_Act(1)) / vP;
    tD_Act(4) = (R_Act(4)-R_Act(1)) / vP;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN SIGNAL PROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Constructing estimated time delays
    %error = tADC;
    error = 2.09964516617352e-07;
    
    for i=2:4;
        tD_Est(i) = tD_Act(i) + (2*rand()-1) * error;
    end

    % Calculating the estimated Time-Of-Arrival
    TOA_Est = (tD_Est(3)^2-tD_Est(2)^2-tD_Est(4)^2) / ...
        (2*(tD_Est(2)-tD_Est(3)+tD_Est(4)));
    
    % Calculating estimated sphere radii
    R_Est(1) = vP*(TOA_Est);
    R_Est(2) = vP*(TOA_Est+tD_Est(2));
    R_Est(3) = vP*(TOA_Est+tD_Est(3));
    R_Est(4) = vP*(TOA_Est+tD_Est(4));
    
    % Determining the estimated source location
    Ping_Est(1) = (R_Est(4)^2-R_Est(1)^2-D^2)/(2*D);
    Ping_Est(2) = (R_Est(2)^2-R_Est(1)^2-D^2)/(2*D);
    Ping_Est(3)  = sqrt(R_Est(1)^2-Ping_Est(1)^2-Ping_Est(2)^2);
    Ping_Est(4)  = -Ping_Est(3);
    
    % Calculating estimated azimuths to source
    if ( isreal(Ping_Est(1)) && isreal(Ping_Est(2)) && isreal(Ping_Est(3)) )
        azimuthH_Est = wrapTo2Pi(atan2(Ping_Est(2),Ping_Est(1))) * (180/pi);
    else
        azimuthH_Est = 0;
        Ping_Est(1) = -1;
        Ping_Est(2) = -1;
        Ping_Est(3) = -1;    
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROCESSING DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    DATA_AZ(1,trialCount) = azimuthH_Act - azimuthH_Est;
end

string1 = sprintf('f_{ADC} = %0.0f [kHz]', fADC/1E+3);
string2 = sprintf('D = %0.1f\\lambda', D/lambda);
string3 = sprintf('Error Analysis for f_{Ping} = %0.1f [kHz]', fPing/1E+3);

figure(1)
        boxplot(DATA_AZ(1,:),'labels',{'Horizontal'});
        hold on;
        line([-1,1],[5,5]);
        line([-1,1],[-5,-5]);
        ylabel('Azimuth Error (deg)');
        legend({string1,string2});
        title(string3);
        hold off;