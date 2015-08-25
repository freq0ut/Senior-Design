%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:  Joshua Simmons
% Started: August 23, 2015
% Status:  Completed
%
% Description:
%
% PINGER ASSUMED TO ALWAYS BE ON!
%
% Uses Time-Difference-of-Arrival (TDOA) to determine the horizontal
% and elevation azimuths to a 30 kHz SINE wave underwater.
%
% 3D Cartesian co-ordinate system. Sensor geometry is a perfect
% tetrahedron.
%
% The origin is best described by a picture:
%   http://i.stack.imgur.com/YAd7z.gif
%
% Coordinates
%   c1 = (d,0,0)            d = D / sqrt(2)
%   c2 = (0,d,0)
%   c3 = (0,0,d)
%   c4 = (d,d,d)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% INITIALIZATION OF PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Global Simulation Parameters
trialTotal = 10; % Total number of iterations of main loop
dwellTime = 0;   % Delay after 1 complete iteration of main loop
fig1_On = false; % Turn on/off visual containing raw time signals and XCs
fig2_On = true;  % Turn on/off visual containing compass and source grid

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
azimuthHs = zeros(1,10);  % Median horizontal azimuth array
azimuthVs = zeros(1,10);  % Median vertical azimuth array
DATA_RAW  = zeros(4,N0);  % Raw data
DATA_CLEAN = zeros(4,N0); % Cleaned data
tD_Act = [0;0;0;0]; % Actual time delays
tD_Est = [0;0;0;0]; % Estimated time delays (Trapezoidal Rule)
XCORR2i = ceil(D/(vP*tADC)); % XCORR2 indices

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
    
    % Time array [s]
    t = 0:tADC:(N0-1)*tADC;
    
    % Calculating random DC offsets
    DC_Offset(1) =  8;
    DC_Offset(2) =  6;
    DC_Offset(3) =  4;
    DC_Offset(4) =  2;
    
    % Constructing the input signals from the time delays and DC offsets
    DATA_RAW(1,:) = DC_Offset(1) + (1.2+0.2*rand())*cos(2*pi*fPing*(t+tD_Act(1))); % Channel 1 (reference)
    DATA_RAW(2,:) = DC_Offset(2) + (1.2+0.2*rand())*cos(2*pi*fPing*(t+tD_Act(2))); % Channel 2
    DATA_RAW(3,:) = DC_Offset(3) + (1.2+0.2*rand())*cos(2*pi*fPing*(t+tD_Act(3))); % Channel 3
    DATA_RAW(4,:) = DC_Offset(4) + (1.2+0.2*rand())*cos(2*pi*fPing*(t+tD_Act(4))); % Channel 4
    
    % Adding TOA to the input signals
    chan = 1;
    while (chan <= 4)
        for i=1:N0;
            if ( i <= round( (TOA_Act+tD_Act(chan))/tADC) );
                if chan == 1
                    DATA_RAW(chan,i) = DC_Offset(1);
                elseif chan == 2
                    DATA_RAW(chan,i) = DC_Offset(2);
                elseif chan == 3
                    DATA_RAW(chan,i) = DC_Offset(3);
                elseif chan == 4
                    DATA_RAW(chan,i) = DC_Offset(4);
                end
            end
        end
        
        chan = chan+1;
    end

    % Adding White Gaussian Noise
    DATA_RAW = awgn(DATA_RAW,SNR);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN SIGNAL PROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Removing DC Offsets    
    chan = 1;
    while (chan <= 4)
        DC_Offset = AVERAGE(t,DATA_RAW(chan,:));
        
        for i=1:N0;
            DATA_CLEAN(chan,i) = DATA_RAW(chan,i) - DC_Offset;        
        end
        
        chan = chan+1;
    end   

    % Determining the estimated time delays using Trapezoidal Rule
    [XC12, XC12_Lags] = XCORR2( DATA_CLEAN(1,:), DATA_CLEAN(2,:), XCORR2i);
    [~,x] = MAXIMUM(XC12_Lags,XC12);
    tD_Est(2) = XC12_Lags(x)*tADC;
    
    [XC13, XC13_Lags] = XCORR2( DATA_CLEAN(1,:), DATA_CLEAN(3,:), XCORR2i);
    [~,x] = MAXIMUM(XC13_Lags,XC13);
    tD_Est(3) = XC13_Lags(x)*tADC;
    
    [XC14, XC14_Lags] = XCORR2( DATA_CLEAN(1,:), DATA_CLEAN(4,:), XCORR2i);
    [~,x] = MAXIMUM(XC14_Lags,XC14);
    tD_Est(4) = XC14_Lags(x)*tADC;
    
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
    
    % Running Medians
    azimuthHs(mod(trialCount,10)+1) = azimuthH_Est;
    azimuthVs(mod(trialCount,10)+1) = azimuthV_Est;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% VISUALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if (fig1_On == true)
        % Raw time signal plots and XCs
        figure(1)
            subplot(2,2,1);
                stem(t*1e6,DATA_RAW(1,:),'-b');
                hold on;
                stem(t*1e6,DATA_RAW(2,:),'-r');
                stem(t*1e6,DATA_RAW(3,:),'-m');
                stem(t*1e6,DATA_RAW(4,:),'-g');
                xlabel('Time [\mus]');
                ylabel('Amplitude');
                legend({'Chan1','Chan2','Chan3','Chan4'});
                string111 = sprintf('f_{samp} = %0.2f [MHz]', fADC/1e6);
                string112 = sprintf('SNR = %0.0f [dB]', SNR);
                title({string111,string112});
                hold off;
            subplot(2,2,2);         
                stem(XC12_Lags*tADC*1e6,XC12,'r');
                hold on;       
                plot(tD_Act(2)*1e6,0,'b.','MarkerSize',20);
                plot(tD_Est(2)*1e6,0,'k.','MarkerSize',20);
                plot(0,0,'w.','MarkerSize',1); 
                string121 = sprintf('td2_{Act} = %f [us]', tD_Act(2)*1e6);
                string122 = sprintf('td2_{Est} = %f [us]', tD_Est(2)*1e6);
                string123 = sprintf('\\Delta td2 = %f [us]', ...
                    (tD_Act(2)-tD_Est(2))*1e6);
                legend({'',string121,string122,string123});
                title('XC_{12}');
                xlabel('Time [\mus]');
                hold off;
            subplot(2,2,3);
                stem(XC13_Lags*tADC*1e6,XC13,'m');
                hold on;
                plot(tD_Act(3)*1e6,0,'b.','MarkerSize',20);
                plot(tD_Est(3)*1e6,0,'k.','MarkerSize',20);
                plot(0,0,'w.','MarkerSize',1);
                string131 = sprintf('td3_{Act} = %f [us]', tD_Act(3)*1e6);
                string132 = sprintf('td3_{Est} = %f [us]', tD_Est(3)*1e6);
                string133 = sprintf('\\Delta td3 = %f [us]', ...
                    (tD_Act(3)-tD_Est(3))*1e6);
                legend({'',string131,string132,string133});
                title('XC_{13}');
                xlabel('Time [\mus]');
                hold off;
            subplot(2,2,4);
                stem(XC14_Lags*tADC*1e6,XC14,'g');
                hold on;
                plot(tD_Act(4)*1e6,0,'b.','MarkerSize',20);
                plot(tD_Est(4)*1e6,0,'k.','MarkerSize',20);
                plot(0,0,'w.','MarkerSize',1);
                string141 = sprintf('td4_{Act} = %f [us]', tD_Act(4)*1e6);
                string142 = sprintf('td4_{Est} = %f [us]', tD_Est(4)*1e6);
                string143 = sprintf('\\Delta td4 = %f [us]', ...
                    (tD_Act(4)-tD_Est(4))*1e6);
                legend({'',string141,string142,string143});
                title('XC_{14}');
                xlabel('Time [\mus]');
                hold off;
    end
        
    if (fig2_On == true)
        
        stringTrials = sprintf('Trial %0.0f / %0.0f', trialCount, trialTotal);
        
        % Compass and source location plots 
        figure(2)
            subplot(2,2,1);
                compass([0 Ping_Act(1)],[0 Ping_Act(2)],'-r');
                hold on;
                scalarXY = sqrt( Ping_Act(1)^2 + Ping_Act(2)^2 ) / ...
                    sqrt( Ping_Est(1)^2 + Ping_Est(2)^2 );
                compass([0 scalarXY*Ping_Est(1)],[0 scalarXY*Ping_Est(2)],'-b');
                %view([90, -90]);
                string211 = sprintf('Actual: %0.1f (deg)', azimuthH_Act);
                string212 = sprintf('Median Est: %0.1f (deg)', median(azimuthHs));
                title({stringTrials,'Horizontal Azimuth',string211,string212,''});
                hold off;
                
            subplot(2,2,3);
                compass([0 Ping_Act(1)],[0 Ping_Act(3)],'-r');
                hold on;
                scalarXZ = sqrt( Ping_Act(1)^2 + Ping_Act(3)^2 ) / ...
                    sqrt( Ping_Est(1)^2 + Ping_Est(3)^2 );
                compass([0 scalarXZ*Ping_Est(1)],[0 scalarXZ*Ping_Est(3)],'-b');
                %view([90, -90]);
                string221 = sprintf('Actual: %0.1f (deg)', azimuthV_Act);
                string222 = sprintf('Median Est: %0.1f (deg)', median(azimuthVs));
                title({stringTrials,'Vertical Azimuth',string221,string222,''});
                hold off;
                
            subplot(2,2,2);
                plot(Ping_Act(1),Ping_Act(2),'r.','MarkerSize',20);
                hold on;
                line([0,Ping_Act(1)],[0,Ping_Act(2)],'Color',[1,0,0]);
                line([0,Ping_Est(1)],[0,Ping_Est(2)],'Color',[0,1,0]);
                grid on;
                xlim([-pingMaxDist,pingMaxDist]);
                ylim([-pingMaxDist,pingMaxDist]);
                hold off;
            subplot(2,2,4);
                plot(Ping_Act(1),Ping_Act(3),'r.','MarkerSize',20);
                hold on;
                line([0,Ping_Act(1)],[0,Ping_Act(3)],'Color',[1,0,0]);
                line([0,Ping_Est(1)],[0,Ping_Est(3)],'Color',[0,1,0]);
                grid on;
                xlim([-pingMaxDist,pingMaxDist]);
                ylim([-pingMaxDist,pingMaxDist]);
                title('XZ Plane');
                hold off;             
    end    
    pause(dwellTime);
end
