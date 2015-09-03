%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:  Joshua Simmons
% Started: August, 2015
% Status:  Incomplete, need to incorporate pinger turning on and off as
%          well as an ADC synchronizing algorithm.
%
% Description:
%
% PINGER ASSUMED TO BE INTERMITTENT AT FIXED INTERVALS!!!
%
% Uses Time-Difference-of-Arrival (TDOA) to determine the azimuth to a
% 30 kHz SINE wave underwater.
%
% 3D Cartesian co-ordinate system with the origin centered in the middle
% of the sensor array. Sensor geometry is square shaped residing all in
% the XY-plane.
%
%   Sensor layout
%
%         (Top View)
%   ---------------------
%   |                   |
%   |         1         |
%   |                   |
%   |     4       2     |
%   |                   |
%   |         3         |
%   |                   |
%   ---------------------
%
% Coordinates
%   N = ( 0, d,0)              D = sqrt(2)*d
%   E = ( d, 0,0)
%   S = ( 0,-d,0)
%   W = (-d, 0,0)
%   P  = (xP,yP,zP)
%
% Sequence of Events
%  1. Initialization of parameters.
%  2. Source location moves in a predictable manner. The actual time delays
%     are computed from the source position.
%  3. Input signals are constructed using the actual time delays. White
%     Gaussian noise is added along with random DC offsets.
%  4. DC Offsets are removed.
%  5. Cross-correlations (XC) are computed for chan2, chan3, and chan4
%     using chan1 as the reference.
%  6. The maximum y-coordinate of each XC is found and the corresponding
%     x-coordinate is multiplied by the sample time. This is the
%     estimated time delay.
%  7. The time delays are plugged into formulas to find the grid
%     coordinates of the source.
%  8. Once the grid coordinates of the source are found, the horizontal and 
%     vertical azimuths are computed.
%  9. Results are visualized.
%
% Be sure that the support functions are in the same directory as this 
% file. Or what you can do is add an extra path to the folder where the
% support functions are located on your PC. You can do this using the 
% "addpath" MatLab command.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;
clc;

addpath('C:\Users\Joshua Simmons\Desktop\Senior_Design\Senior-Design\MATLAB\Support_Functions');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETER INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Global Simulation Parameters
trialTotal = 1E+3;    % Total number of iterations of main loop
dwellTime = 0;        % Delay after 1 complete iteration of main loop
fig2_On = true;       % Turn on/off visual containing compass and source grid
OS_dwellTime = false; % Gives you time to make to go full screen at the start
                      % of the simulation haha!

% Pinger Properties
SNR  = 20;         % Signal to Noise Ratio [dB]
fPing = 30E+3;     % Source freq [Hz]
tPing = 1/fPing;   % Source period [s]
vP = 1482;         % Propagation Velocity [m/s]
lambda = vP/fPing; % Wavelength [m]
pingMaxDist = 0.5; % Pinger max distance from sensors [m]

% Hydrophone Properties
D = lambda;      % Hydrophone spacing [m]
d = D / sqrt(2); % For coordinates of sensors

% ADC
fADC = 1800E+3; % Sample freq [Hz]
tADC = 1/fADC;  % Sample period [s]
N0 = 2^11;      % Samples per frame

% Microcontroller Properties
azimuthH_Est  = 0; % For 1st iteration
azimuthV_Est  = 0; % For 1st iteration
azimuthV2_Est = 0; % For 1st iteration
azimuthHs =  zeros(1,10); % Median horizontal azimuth array
azimuthVs =  zeros(1,10); % Median vertical azimuth array
azimuthV2s =  zeros(1,10); % Second Median vertical azimuth array
DATA = zeros(4,N0);  % Raw data
tD_Act  = [0;0;0;0]; % Actual time delays
tD_Est = [0;0;0;0];  % Best guess as to estimated time delays
tD_EstP = [0;0;0;0]; % Primary estimated time delays
tD_EstS =  zeros(4,ceil(2*D/lambda)); % Secondary estimated time delays
TOA_Est = zeros(1,4); % Estimated Time-Of-Arrivals
XCORR2i = ceil(sqrt(2)*D/(vP*tADC)); % XCORR2 indices
MPD = 60;  % Minimum Peak Distance
xNBRS = 2; % Neighbors to look to the left and right of
THD = 0.8; % Threshold

% Enabling GPU Computing
%DATA = gpuArray(DATA);

% Testing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTRUCTING INPUT SIGNALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% START MAIN LOOP
for trialCount = 1:trialTotal;
    
    % Actual pinger location (xP,yP,zP)
    if ( mod(trialCount,20) == 1 )
        Ping_Act(1) =  (2*rand()-1)*(pingMaxDist-D) + D;
        Ping_Act(2) =  (2*rand()-1)*(pingMaxDist-D) + D;
        Ping_Act(3) =  (2*rand()-1)*(pingMaxDist-D) + D;
    end
    
    % Actual azimuth to pinger
    azimuthH_Act = wrapTo2Pi(atan2(Ping_Act(2),Ping_Act(1))) * (180/pi);
    azimuthV_Act = wrapTo2Pi(atan2(Ping_Act(3),Ping_Act(1))) * (180/pi);
    
    % Actual sphere radii
    R_Act(1) = sqrt( (Ping_Act(1)  )^2 + (Ping_Act(2)-d)^2 + (Ping_Act(3)  )^2 );
    R_Act(2) = sqrt( (Ping_Act(1)-d)^2 + (Ping_Act(2)  )^2 + (Ping_Act(3)  )^2 );
    R_Act(3) = sqrt( (Ping_Act(1)  )^2 + (Ping_Act(2)+d)^2 + (Ping_Act(3)  )^2 );
    R_Act(4) = sqrt( (Ping_Act(1)+d)^2 + (Ping_Act(2)  )^2 + (Ping_Act(3)  )^2 );
    
    % Actual Time-Of-Arrival
    TOA_Act = R_Act(1)/vP;

    % Actual time delays
    tD_Act(2) = (R_Act(2)-R_Act(1)) / vP;
    tD_Act(3) = (R_Act(3)-R_Act(1)) / vP;
    tD_Act(4) = (R_Act(4)-R_Act(1)) / vP;
    
    % Time array [s]
    t = 0:tADC:(N0-1)*tADC;
    
    % DC Offsets
    DC_Offset(1) =  8;
    DC_Offset(2) =  6;
    DC_Offset(3) =  4;
    DC_Offset(4) =  2;
    
    % Incorporating DC offsets and time delays
    DATA(1,:) = DC_Offset(1) + (1.2+0.2*rand())*cos(2*pi*fPing*(t+tD_Act(1))); % Channel 1
    DATA(2,:) = DC_Offset(2) + (1.2+0.2*rand())*cos(2*pi*fPing*(t+tD_Act(2))); % Channel 2
    DATA(3,:) = DC_Offset(3) + (1.2+0.2*rand())*cos(2*pi*fPing*(t+tD_Act(3))); % Channel 3
    DATA(4,:) = DC_Offset(4) + (1.2+0.2*rand())*cos(2*pi*fPing*(t+tD_Act(4))); % Channel 4
    
    % Incorporating TOA Actual
    for i=1:round( (TOA_Act+tD_Act(1))/tADC );
        DATA(1,i) = DC_Offset(1);
    end

    for i=1:round( (TOA_Act+tD_Act(2))/tADC );
        DATA(2,i) = DC_Offset(2);
    end

    for i=1:round( (TOA_Act+tD_Act(3))/tADC );
        DATA(3,i) = DC_Offset(3);
    end

    for i=1:round( (TOA_Act+tD_Act(4))/tADC );
        DATA(4,i) = DC_Offset(4);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN SIGNAL PROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    chan = 1;
    while (chan <= 4)
        % Removing DC offsets
        DC_Offset = AVERAGE(t,DATA(chan,:));
        
        for i=1:N0;
            DATA(chan,i) = DATA(chan,i) - DC_Offset;        
        end
        
        % Estimated TOAs
        iBreak = BREAK_THRESHOLD(DATA(chan,:),THD);
        TOA_Est(chan) = iBreak*tADC;
        
        chan = chan+1;
    end
    
    chan = 1;
    while (chan <= 4)
        % Primary estimated time delays
        tD_EstP(chan) = TOA_Est(chan) - TOA_Est(1);
        
        % Secondary estimated time delays
        [XC, XC_Lags] = XCORR3( DATA(1,:), DATA(chan,:), XCORR2i );
        [~,pkLocs] = FIND_PEAKS(XC,MPD,xNBRS);

        for i=1:length(pkLocs);
            tD_EstS(chan,i) = XC_Lags(pkLocs(i))*tADC;
        end
        
        % Esimated time delays
        tD_Est(chan) = COMPARE( tD_EstP(chan), tD_EstS(chan,:) );
        
        chan = chan + 1;
    end

    % Estimated Time-Of-Arrival
    TOA_Est_Scalar = (tD_Est(2)^2+tD_Est(4)^2-tD_Est(3)^2) / ...
        (2*(tD_Est(3)-tD_Est(2)-tD_Est(4)));
    
    % Estimated sphere radii
    R_Est(1) = vP*(TOA_Est_Scalar);
    R_Est(2) = vP*(TOA_Est_Scalar+tD_Est(2));
    R_Est(3) = vP*(TOA_Est_Scalar+tD_Est(3));
    R_Est(4) = vP*(TOA_Est_Scalar+tD_Est(4));
    
    % Estimated pinger location (xP,yP,zP)
    Ping_Est(1) = (R_Est(4)^2-R_Est(2)^2)/(4*d);
    Ping_Est(2) = (R_Est(3)^2-R_Est(1)^2)/(4*d);
    Ping_Est(3)  = sqrt(R_Est(1)^2-Ping_Est(1)^2-(Ping_Est(2)-d)^2);
    Ping_Est(4)  = -Ping_Est(3);
    
    % Estimated azimuths to pinger
    if ( isreal(Ping_Est(1)) && isreal(Ping_Est(2)) && isreal(Ping_Est(3)) ...
            && ~isnan(Ping_Est(1)) && ~isnan(Ping_Est(2)) && ~isnan(Ping_Est(3)))

        azimuthH_Est  = wrapTo2Pi(atan2( Ping_Est(2),Ping_Est(1))) * (180/pi);
        azimuthV_Est  = wrapTo2Pi(atan2( Ping_Est(3),Ping_Est(1))) * (180/pi);
        azimuthV2_Est = wrapTo2Pi(atan2(-Ping_Est(3),Ping_Est(1))) * (180/pi);
        
        % Running medians get updated with new information
        azimuthHs (mod(trialCount,10)+1) = azimuthH_Est;
        azimuthVs (mod(trialCount,10)+1) = azimuthV_Est;
        azimuthV2s(mod(trialCount,10)+1) = azimuthV2_Est;
    else
        % Running medians get updated with old information
        azimuthHs (mod(trialCount,10)+1) = azimuthH_Est;
        azimuthVs (mod(trialCount,10)+1) = azimuthV_Est;
        azimuthV2s(mod(trialCount,10)+1) = azimuthV2_Est;
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VISUALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if (fig2_On == true)
        
        stringTrials = sprintf('Trial %0.0f / %0.0f', trialCount, trialTotal);
        
        % Compass and source location plots 
        figure(2)
            subplot(2,2,1);
                compass([0 Ping_Act(1)],[0 Ping_Act(2)],'-k');
                hold on;
                scalarXY = sqrt( Ping_Act(1)^2 + Ping_Act(2)^2 ) / ...
                    sqrt( Ping_Est(1)^2 + Ping_Est(2)^2 );
                compass([0 scalarXY*Ping_Est(1)],[0 scalarXY*Ping_Est(2)],'-c');
                %view([90, -90]);
                string211 = sprintf('Actual: %0.1f (deg)', azimuthH_Act);
                string212 = sprintf('Median Est: %0.1f (deg)', median(azimuthHs));
                title({stringTrials,'Horizontal Azimuth',string211,string212,''});
                hold off;
                
            subplot(2,2,3);
                compass([0 Ping_Act(1)],[0 Ping_Act(3)],'-k');
                hold on;
                scalarXZ = sqrt( Ping_Act(1)^2 + Ping_Act(3)^2 ) / ...
                    sqrt( Ping_Est(1)^2 + Ping_Est(3)^2 );
                compass([0 scalarXZ*Ping_Est(1)],[0  scalarXZ*Ping_Est(3)],'-c');
                compass([0 scalarXZ*Ping_Est(1)],[0 -scalarXZ*Ping_Est(3)],'-c');
                %view([90, -90]);
                string221 = sprintf('Actual: %0.1f (deg)', azimuthV_Act);
                string222 = sprintf('Median Est: %0.1f (deg)', median(azimuthVs));
                string223 = sprintf('Median Est: %0.1f (deg)', median(azimuthV2s));
                title({stringTrials,'Vertical Azimuth',string221,string222,string223,''});
                hold off;
                
            subplot(2,2,2);
                plot( 0, d,'b.','MarkerSize',20);
                hold on;
                plot( d, 0,'r.','MarkerSize',20);
                plot( 0,-d,'m.','MarkerSize',20);
                plot(-d, 0,'g.','MarkerSize',20);
                plot(Ping_Act(1),Ping_Act(2),'k.','MarkerSize',20);
                plot(Ping_Est(1),Ping_Est(2),'c.','MarkerSize',20);
                %line([0,Ping_Act(1)],[0,Ping_Act(2)],'Color',[1,0,0]);
                ezpolar(@(x)N0*tADC*vP);
                grid on;
                xlim([-2*pingMaxDist,2*pingMaxDist]);
                ylim([-2*pingMaxDist,2*pingMaxDist]);
                legend({'chan1','chan2','chan3','chan4'});
                title('XY Plane');
                hold off;
            subplot(2,2,4);
                plot(Ping_Act(1),Ping_Act(3),'k.','MarkerSize',20);
                hold on;
                plot(Ping_Est(1), Ping_Est(3),'c.','MarkerSize',20);
                plot(Ping_Est(1),-Ping_Est(3),'c.','MarkerSize',20);
                %line([0,Ping_Act(1)],[0,Ping_Act(3)],'Color',[1,0,0]);
                ezpolar(@(x)N0*tADC*vP);
                grid on;
                xlim([-2*pingMaxDist,2*pingMaxDist]);
                ylim([-2*pingMaxDist,2*pingMaxDist]);
                title('XZ Plane');
                hold off; 
    end
    
    if (OS_dwellTime == false)
        pause(10);
        OS_dwellTime = true;
    end
    
    pause(dwellTime);
    
end
