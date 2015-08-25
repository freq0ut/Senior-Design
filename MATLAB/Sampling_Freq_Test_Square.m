close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% INITIALIZATION OF PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Global Simulation Parameters
trialTotal = 1e4;  % Total number of iterations of main loop

% ADC
fS = 1800e3;
tS = 1/fS;

% Source Properties
SNR  = 10;        % Signal to Noise Ratio [dB]
fSce = 30e3;      % Source freq [Hz]
Tsce = 1/fSce;    % Source period [s]
vP   = 1482;      % Propagation Velocity [m/s]
lambda = vP/fSce; % Wavelength [m]
S_Act = [0;0;0];  % Initialization of source location [x,y,z] in [m]

% Hydrophone Properties
D = lambda;       % Hydrophone spacing [m]

% Azimuths Data
DATA_Azimuths = zeros(2,trialTotal);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% CONTRUCTING INPUT SIGNALS %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% START MAIN LOOP
for trialCount = 1:trialTotal;
    % Simulating a source moving in the water.
    S_Act(1) = 100*(2*rand()-1);
    S_Act(2) = 100*(2*rand()-1);
    S_Act(3) = 100*(2*rand()-1);

    % Calculating actual azimuths to source
    azimuthH_Act = wrapTo2Pi(atan2(S_Act(2),S_Act(1))) * (180/pi);
    azimuthV_Act = wrapTo2Pi(atan2(S_Act(3),S_Act(1))) * (180/pi);
    
    % Determining the actual time delays
    R1_Act = sqrt( (S_Act(1)  )^2 + (S_Act(2)  )^2 + (S_Act(3)  )^2 );
    R2_Act = sqrt( (S_Act(1)  )^2 + (S_Act(2)+D)^2 + (S_Act(3)  )^2 );
    R3_Act = sqrt( (S_Act(1)+D)^2 + (S_Act(2)+D)^2 + (S_Act(3)  )^2 );
    R4_Act = sqrt( (S_Act(1)+D)^2 + (S_Act(2)  )^2 + (S_Act(3)  )^2 );
    
    tD_Act = [0;0;0;0];
    TOA_Act = R1_Act/vP;
    tD_Act(2) = (R2_Act-R1_Act) / vP;
    tD_Act(3) = (R3_Act-R1_Act) / vP;
    tD_Act(4) = (R4_Act-R1_Act) / vP;
    
    % Constructing estimated time delays
    error = 1/(2*10e6);
    tD_Est = tD_Act;
    
    for i=2:4;
        tD_Est(i) = tD_Est(i) + (2*rand()-1) * error;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN SIGNAL PROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    S_Est = [0;0;0;0];

    % Calculating the estimated Time-Of-Arrival
    TOA_Est = (tD_Est(3)^2-tD_Est(2)^2-tD_Est(4)^2) / ...
        (2*(tD_Est(2)-tD_Est(3)+tD_Est(4)));
    
    % Calculating sphere radii
    R1_Est = vP*(TOA_Est);
    R2_Est = vP*(TOA_Est+tD_Est(2));
    R3_Est = vP*(TOA_Est+tD_Est(3));
    R4_Est = vP*(TOA_Est+tD_Est(4));
    
    % Determining source location
    S_Est(1) = (R4_Est^2-R1_Est^2-D^2)/(2*D);
    S_Est(2) = (R2_Est^2-R1_Est^2-D^2)/(2*D);
    S_Est(3)  = sqrt(R1_Est^2-S_Est(1)^2-S_Est(2)^2);
    S_Est(4)  = -S_Est(3);
    
    % Calculating estimated azimuths to source
    if ( isreal(S_Est(1)) && isreal(S_Est(2)) && isreal(S_Est(3)) )
        azimuthH_Est = wrapTo2Pi(atan2(S_Est(2),S_Est(1))) * (180/pi);
        azimuthV_Est = wrapTo2Pi(atan2(S_Est(3),S_Est(1))) * (180/pi);
        
        DATA_Azimuths(1,trialCount) = azimuthH_Act - azimuthH_Est;
        DATA_Azimuths(2,trialCount) = azimuthV_Act - azimuthV_Est;
    else
        trialCount = trialCount - 1;
    end
end

figure(1)
    boxplot(DATA_Azimuths(1,:),'labels',{'H'});
    hold on;
    line([-1,1],[5,5]);
    line([-1,1],[-5,-5]);
    hold off;
    
% 5 degrees
% D = 2*lambda, fS =   1.8e6 EXCELLENT
% D = lambda/3, fS = > 2 GHz