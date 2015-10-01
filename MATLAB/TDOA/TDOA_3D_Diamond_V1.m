%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:  Joshua Simmons
% Started: August, 2015
% Status:  COMPLETED
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

%addpath('C:\Users\Joshua Simmons\Desktop\Senior_Design\Senior-Design\MATLAB\Support_Functions');
addpath('/Users/betio32/Documents/myGitHub/Senior-Design/MATLAB/Support_Functions');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETER INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Global Simulation Parameters
trialTotal = 1;    % Total number of iterations of main loop
dwellTime = 0.5;      % Delay after 1 complete iteration of main loop
fig1_On = false;      % Turn on/off visual containing raw time signals and XCs
fig2_On = true;       % Turn on/off visual containing compass and source grid
fig3_On = false;      % Turn on/off visual containing RAW FFT and iFFT
fig4_On = false;      % Turn on/off visual containing CLEAN FFT and iFFT
OS_dwellTime = false; % Gives you time to make to go full screen at the start
                      % of the simulation haha!

% Pinger Properties
SNR  = 0;         % Signal to Noise Ratio [dB]
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
N0 = 2^10;      % Samples per frame

% FFT
t0 = N0*tADC;   % Truncation Time Interval [s]
f0 = 1/t0;      % Frequency Resolution [Hz]
f = f0*(-N0/2:N0/2-1); % Double Sided Freq Array [Hz]

% Digital Ideal BPF
H = zeros(1,N0);
fCent = 30E+3;
chanHalfBW = 5E+3;
for i = 1:N0; 
    if ( abs(f(i)) >= fCent-chanHalfBW && abs(f(i)) <= fCent+chanHalfBW )
        H(i) = 1;
    end
end

% Microcontroller Properties
azimuthH_Est  = 0; % For 1st iteration
azimuthV_Est  = 0; % For 1st iteration
azimuthV2_Est = 0; % For 1st iteration
azimuthHs = zeros(1,10); % Median horizontal azimuth array
azimuthVs = zeros(1,10); % Median vertical azimuth array
azimuthV2s = zeros(1,10); % Second Median vertical azimuth array
DATA_RAW_t = zeros(4,N0); % Raw time data
DATA_RAW_f = zeros(4,N0); % Raw frequency data
DATA_CLEAN_t = zeros(4,N0); % Cleaned time data
DATA_CLEAN_f = zeros(4,N0); % Cleaned frequency data
tD_Act  = [0;0;0;0]; % Actual time delays
tD_EstP = [0;0;0;0]; % Primary estimated time delays
tD_EstS =  zeros(4,ceil(2*D/lambda)); % Secondary estimated time delays
TOA_Est = zeros(1,4); % Estimated Time-Of-Arrivals
XCORR2i = ceil(sqrt(2)*D/(vP*tADC)); % XCORR2 indices
MPD = 60;  % Minimum Peak Distance
xNBRS = 2; % Neighbors to look to the left and right of
THD = 0.6; % Threshold

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTRUCTING INPUT SIGNALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% START MAIN LOOP
for trialCount = 1:trialTotal;
    
    % Actual pinger location (xP,yP,zP)
    if ( mod(trialCount,20) == 1) 
        Ping_Act(1) =  (2*rand()-1)*(pingMaxDist-D) + D;
        Ping_Act(2) =  (2*rand()-1)*(pingMaxDist-D) + D;
        Ping_Act(3) =  (2*rand()-1)*(pingMaxDist-D) + D;
    end
    
    % Actual azimuth to pinger
    azimuthH_Act = WRAPTO2PI(atan2(Ping_Act(2),Ping_Act(1))) * (180/pi);
    azimuthV_Act = WRAPTO2PI(atan2(Ping_Act(3),Ping_Act(1))) * (180/pi);
    
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
    
    % Incorporating time delays
    DATA_RAW_t(1,:) = cos(2*pi*fPing*(t+tD_Act(1))); % Channel 1
    DATA_RAW_t(2,:) = cos(2*pi*fPing*(t+tD_Act(2))); % Channel 2
    DATA_RAW_t(3,:) = cos(2*pi*fPing*(t+tD_Act(3))); % Channel 3
    DATA_RAW_t(4,:) = cos(2*pi*fPing*(t+tD_Act(4))); % Channel 4
            
    % Incorporating TOA Actual (POST)
    for i=1:round( (TOA_Act+tD_Act(1))/tADC );
        DATA_RAW_t(1,i) = 0;
    end

    for i=1:round( (TOA_Act+tD_Act(2))/tADC );
        DATA_RAW_t(2,i) = 0;
    end

    for i=1:round( (TOA_Act+tD_Act(3))/tADC );
        DATA_RAW_t(3,i) = 0;
    end

    for i=1:round( (TOA_Act+tD_Act(4))/tADC );
        DATA_RAW_t(4,i) = 0;
    end

%     % Fast-Fourier Transform (Cooley-Tukey) into Frequency Domain
%     DATA_RAW_f(1,:) = fftshift(fft(DATA_RAW_t(1,:))) / N0;
%     DATA_RAW_f(2,:) = fftshift(fft(DATA_RAW_t(2,:))) / N0;
%     DATA_RAW_f(3,:) = fftshift(fft(DATA_RAW_t(3,:))) / N0;
%     DATA_RAW_f(4,:) = fftshift(fft(DATA_RAW_t(4,:))) / N0;
%     
%     % Adding low frequency ( <= 10 kHz) contamination
%     i = 1;
%     while (i<=N0)
%         if ( f(i) >= -10E+3 )
%            iLow = i;
%            i = N0;
%         end
%         
%         i = i+1;
%     end
%     
%     i = N0;
%     while (i>=1)
%         if ( f(i) <= 10E+3 )
%            iHigh = i;
%            i = 1;
%         end
%         
%         i = i-1;
%     end
%     
%     DATA_RAW_f(1,iLow:iHigh) = awgn( DATA_RAW_f(1,iLow:iHigh),0);
%     DATA_RAW_f(2,iLow:iHigh) = awgn( DATA_RAW_f(2,iLow:iHigh),0);
%     DATA_RAW_f(3,iLow:iHigh) = awgn( DATA_RAW_f(3,iLow:iHigh),0);
%     DATA_RAW_f(4,iLow:iHigh) = awgn( DATA_RAW_f(4,iLow:iHigh),0);
%     
%     % Fast-Fourier Transform (Cooley-Tukey) back into Time Domain
%     DATA_RAW_t(1,:) = ifft(ifftshift(DATA_RAW_f(1,:))) * N0;
%     DATA_RAW_t(2,:) = ifft(ifftshift(DATA_RAW_f(2,:))) * N0;
%     DATA_RAW_t(3,:) = ifft(ifftshift(DATA_RAW_f(3,:))) * N0;
%     DATA_RAW_t(4,:) = ifft(ifftshift(DATA_RAW_f(4,:))) * N0;
    
    % Adding White Gaussian Noise
    DATA_RAW_t = awgn(DATA_RAW_t,SNR);
    
    % Exporting to CSV for simulations
    csvwrite('/Users/betio32/Desktop/TDOA_Chan1_SimData001.CSV', DATA_RAW_t(1,:)');
    csvwrite('/Users/betio32/Desktop/TDOA_Chan2_SimData001.CSV', DATA_RAW_t(2,:)');
    csvwrite('/Users/betio32/Desktop/TDOA_Chan3_SimData001.CSV', DATA_RAW_t(3,:)');
    csvwrite('/Users/betio32/Desktop/TDOA_Chan4_SimData001.CSV', DATA_RAW_t(4,:)');
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN SIGNAL PROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Fast-Fourier Transform (Cooley-Tukey) into Frequency Domain
    DATA_CLEAN_f(1,:) = fftshift(fft(DATA_RAW_t(1,:))) / N0;
    DATA_CLEAN_f(2,:) = fftshift(fft(DATA_RAW_t(2,:))) / N0;
    DATA_CLEAN_f(3,:) = fftshift(fft(DATA_RAW_t(3,:))) / N0;
    DATA_CLEAN_f(4,:) = fftshift(fft(DATA_RAW_t(4,:))) / N0;
    
    % Applying Digital Ideal BPF
    DATA_CLEAN_f(1,:) = H.*DATA_CLEAN_f(1,:);
    DATA_CLEAN_f(2,:) = H.*DATA_CLEAN_f(2,:);
    DATA_CLEAN_f(3,:) = H.*DATA_CLEAN_f(3,:);
    DATA_CLEAN_f(4,:) = H.*DATA_CLEAN_f(4,:);
    
    % Fast-Fourier Transform (Cooley-Tukey) back into Time Domain
    DATA_CLEAN_t(1,:) = ifft(ifftshift(DATA_CLEAN_f(1,:))) * N0;
    DATA_CLEAN_t(2,:) = ifft(ifftshift(DATA_CLEAN_f(2,:))) * N0;
    DATA_CLEAN_t(3,:) = ifft(ifftshift(DATA_CLEAN_f(3,:))) * N0;
    DATA_CLEAN_t(4,:) = ifft(ifftshift(DATA_CLEAN_f(4,:))) * N0;

    % Estimated TOAs
    % LR for heads triggered
    % RL for tails triggered
    iBreak1 = BREAK_THRESHOLD(DATA_CLEAN_t(1,:),THD,'LR');
    TOA_Est(1) = (N0-iBreak1)*tADC;

    iBreak2 = BREAK_THRESHOLD(DATA_CLEAN_t(2,:),THD,'LR');
    TOA_Est(2) = (N0-iBreak2)*tADC;

    iBreak3 = BREAK_THRESHOLD(DATA_CLEAN_t(3,:),THD,'LR');
    TOA_Est(3) = (N0-iBreak3)*tADC;

    iBreak4 = BREAK_THRESHOLD(DATA_CLEAN_t(4,:),THD,'LR');
    TOA_Est(4) = (N0-iBreak4)*tADC;

    % Primary estimated time delays
    tD_EstP(2) = TOA_Est(1) - TOA_Est(2);
    tD_EstP(3) = TOA_Est(1) - TOA_Est(3);
    tD_EstP(4) = TOA_Est(1) - TOA_Est(4);
    
    % Secondary estimated time delays
    [XC12, XC12_Lags] = XCORR2( DATA_CLEAN_t(1,:), DATA_CLEAN_t(2,:), XCORR2i );
    [~,pkLocs12] = FIND_PEAKS(XC12,THD,MPD,xNBRS);

    for i=1:length(pkLocs12);
        tD_EstS(2,i) = XC12_Lags(pkLocs12(i))*tADC;
    end

    [XC13, XC13_Lags] = XCORR2( DATA_CLEAN_t(1,:), DATA_CLEAN_t(3,:), XCORR2i );
    [~,pkLocs13] = FIND_PEAKS(XC13,THD,MPD,xNBRS);

    for i=1:length(pkLocs13);
        tD_EstS(3,i) = XC13_Lags(pkLocs13(i))*tADC;
    end
    
    [XC14, XC14_Lags] = XCORR2( DATA_CLEAN_t(1,:), DATA_CLEAN_t(4,:), XCORR2i );
    [~,pkLocs14] = FIND_PEAKS(XC14,THD,MPD,xNBRS);

    for i=1:length(pkLocs14);
        tD_EstS(4,i) = XC14_Lags(pkLocs14(i))*tADC;
    end
    
    % Estimated time delays
    tD_Est2 = COMPARE( tD_EstP(2), tD_EstS(2,:) );
    tD_Est3 = COMPARE( tD_EstP(3), tD_EstS(3,:) );
    tD_Est4 = COMPARE( tD_EstP(4), tD_EstS(4,:) );

    % Estimated Time-Of-Arrival
    TOA_Est_Scalar = (tD_Est2^2+tD_Est4^2-tD_Est3^2) / ...
        (2*(tD_Est3-tD_Est2-tD_Est4));
    
    % Estimated sphere radii
    R_Est(1) = vP*(TOA_Est_Scalar);
    R_Est(2) = vP*(TOA_Est_Scalar+tD_Est2);
    R_Est(3) = vP*(TOA_Est_Scalar+tD_Est3);
    R_Est(4) = vP*(TOA_Est_Scalar+tD_Est4);
    
    % Estimated pinger location (xP,yP,zP)
    Ping_Est(1) = (R_Est(4)^2-R_Est(2)^2)/(4*d);
    Ping_Est(2) = (R_Est(3)^2-R_Est(1)^2)/(4*d);
    Ping_Est(3)  = sqrt(R_Est(1)^2-Ping_Est(1)^2-(Ping_Est(2)-d)^2);
    Ping_Est(4)  = -Ping_Est(3);
    
    % Estimated azimuths to pinger
    if ( isreal(Ping_Est(1)) && isreal(Ping_Est(2)) && isreal(Ping_Est(3)) ...
            && ~isnan(Ping_Est(1)) && ~isnan(Ping_Est(2)) && ~isnan(Ping_Est(3)))

        azimuthH_Est  = WRAPTO2PI(atan2( Ping_Est(2),Ping_Est(1))) * (180/pi);
        azimuthV_Est  = WRAPTO2PI(atan2( Ping_Est(3),Ping_Est(1))) * (180/pi);
        azimuthV2_Est = WRAPTO2PI(atan2(-Ping_Est(3),Ping_Est(1))) * (180/pi);
        
        % Running medians get updated with new information
        azimuthHs(mod(trialCount,10)+1)  = azimuthH_Est;
        azimuthVs(mod(trialCount,10)+1)  = azimuthV_Est;
        azimuthV2s(mod(trialCount,10)+1) = azimuthV2_Est;
    else
        % Running medians get updated with old information
        azimuthHs(mod(trialCount,10)+1)  = azimuthH_Est;
        azimuthVs(mod(trialCount,10)+1)  = azimuthV_Est;
        azimuthV2s(mod(trialCount,10)+1) = azimuthV2_Est;
    end
    
    % Taking the medians
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VISUALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if (fig1_On == true)
        % Raw time signal plots and XCs
        figure(1)
            subplot(2,2,1);
                plot(t*1E+6,DATA_RAW_t(1,:),'-b');
                hold on;
                plot(t*1E+6,DATA_RAW_t(2,:),'-r');
                plot(t*1E+6,DATA_RAW_t(3,:),'-m');
                plot(t*1E+6,DATA_RAW_t(4,:),'-g');
                xlabel('Time [\mus]');
                ylabel('Amplitude');
                legend({'Chan1','Chan2','Chan3','Chan4'});
                string111 = sprintf('f_{samp} = %0.2f [MHz]', fADC/1E+6);
                string112 = sprintf('SNR = %0.0f [dB]', SNR);
                title({string111,string112});
                hold off;
            subplot(2,2,2);         
                stem(XC12_Lags*tADC*1E+6,XC12,'r');
                hold on;       
                plot(tD_Act(2)*1E+6,0,'k.','MarkerSize',20);
                plot(tD_Est2*1E+6,0,'b.','MarkerSize',20);
                plot(0,0,'w.','MarkerSize',1); 
                string121 = sprintf('td2_{Act} = %f [\\mus]', tD_Act(2)*1E+6);
                string122 = sprintf('td2_{Est} = %f [\\mus]', tD_Est2*1E+6);
                string123 = sprintf('\\Delta td2 = %f [\\mus]', ...
                    (tD_Act(2)-tD_Est2)*1e6);
                legend({'',string121,string122,string123});
                title('XC_{12}');
                xlabel('Time [\mus]');
                hold off;
            subplot(2,2,3);
                stem(XC13_Lags*tADC*1E+6,XC13,'m');
                hold on;
                plot(tD_Act(3)*1E+6,0,'k.','MarkerSize',20);
                plot(tD_Est3*1E+6,0,'b.','MarkerSize',20);
                plot(0,0,'w.','MarkerSize',1);
                string131 = sprintf('td3_{Act} = %f [\\mus]', tD_Act(3)*1E+6);
                string132 = sprintf('td3_{Est} = %f [\\mus]', tD_Est3*1E+6);
                string133 = sprintf('\\Delta td3 = %f [\\mus]', ...
                    (tD_Act(3)-tD_Est3)*1E+6);
                legend({'',string131,string132,string133});
                title('XC_{13}');
                xlabel('Time [\mus]');
                hold off;
            subplot(2,2,4);
                stem(XC14_Lags*tADC*1E+6,XC14,'g');
                hold on;
                plot(tD_Act(4)*1E+6,0,'k.','MarkerSize',20);
                plot(tD_Est4*1E+6,0,'b.','MarkerSize',20);
                plot(0,0,'w.','MarkerSize',1);
                string141 = sprintf('td4_{Act} = %f [\\mus]', tD_Act(4)*1E+6);
                string142 = sprintf('td4_{Est} = %f [\\mus]', tD_Est4*1E+6);
                string143 = sprintf('\\Delta td4 = %f [\\mus]', ...
                    (tD_Act(4)-tD_Est4)*1E+6);
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
                line([-pingMaxDist,pingMaxDist],[pingMaxDist,-pingMaxDist],'Color',[1,0,0]);
                line([-pingMaxDist,pingMaxDist],[-pingMaxDist,pingMaxDist],'Color',[1,0,0]);
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
    
    if (fig3_On == true)
        % RAW FFT and iFFT
        figure(3);
            subplot(2,4,1);
                plot(t*1E+6,DATA_RAW_t(1,:),'-b');
                hold on;
                grid on;
                grid minor;
                xlabel('Time [\mus]');
                ylabel('Amplitude');
                ylim([-10,10]);
                legend('chan1');
                hold off;
            subplot(2,4,2);
                plot(t*1E+6,DATA_RAW_t(2,:),'-r');
                hold on;
                grid on;
                grid minor;
                xlabel('Time [\mus]');
                ylabel('Amplitude');
                ylim([-10,10]);
                legend('chan2');
                hold off;
            subplot(2,4,5);
                plot(t*1E+6,DATA_RAW_t(3,:),'-m');
                hold on;
                grid on;
                grid minor;
                xlabel('Time [\mus]');
                ylabel('Amplitude');
                ylim([-10,10]);
                legend('chan3');
                hold off;
            subplot(2,4,6);
                plot(t*1E+6,DATA_RAW_t(4,:),'-g');
                hold on;
                grid on;
                grid minor;
                xlabel('Time [\mus]');
                ylabel('Amplitude');
                ylim([-10,10]);
                legend('chan4');
                hold off;
        
            subplot(2,4,3);
                stem(f/1E+3,abs(DATA_RAW_f(1,:)),'-b');
                hold on;
                stem(f/1E+3,H,'-y');
                grid on;
                grid minor;
                xlim([2*(-fCent-chanHalfBW)/1E+3,2*(fCent+chanHalfBW)/1E+3]);
                xlabel('Frequency [kHz]');
                ylabel('Magnitude');
                ylim([0,3]);
                legend('chan1');
                hold off;
            subplot(2,4,4);
                stem(f/1E+3,abs(DATA_RAW_f(2,:)),'-r');
                hold on;
                stem(f/1E+3,H,'-y');
                grid on;
                grid minor;
                xlim([2*(-fCent-chanHalfBW)/1E+3,2*(fCent+chanHalfBW)/1E+3]);
                xlabel('Frequency [kHz]');
                ylabel('Magnitude');
                ylim([0,3]);
                legend('chan2');
                hold off;
            subplot(2,4,7);
                stem(f/1E+3,abs(DATA_RAW_f(3,:)),'-m');
                hold on;
                stem(f/1E+3,H,'-y');
                grid on;
                grid minor;
                xlim([2*(-fCent-chanHalfBW)/1E+3,2*(fCent+chanHalfBW)/1E+3]);
                xlabel('Frequency [kHz]');
                ylabel('Magnitude');
                ylim([0,3]);
                legend('chan3');
                hold off;
            subplot(2,4,8);
                stem(f/1E+3,abs(DATA_RAW_f(4,:)),'-g');
                hold on;
                stem(f/1E+3,H,'-y');
                grid on;
                grid minor;
                xlim([2*(-fCent-chanHalfBW)/1E+3,2*(fCent+chanHalfBW)/1E+3]);
                xlabel('Frequency [kHz]');
                ylabel('Magnitude');
                ylim([0,3]);
                legend('chan4');
                hold off;
                
            
    end
    
    if (fig4_On == true)
        % CLEAN FFT and iFFT
        figure(4);
            subplot(2,4,1);
                stem(f/1E+3,abs(DATA_CLEAN_f(1,:)),'-b');
                hold on;
                %stem(f/1E+3,H,'-y');
                grid on;
                grid minor;
                xlim([2*(-fCent-chanHalfBW)/1E+3,2*(fCent+chanHalfBW)/1E+3]);
                xlabel('Frequency [kHz]');
                ylabel('Magnitude');
                legend('chan1');
                hold off;
            subplot(2,4,2);
                stem(f/1E+3,abs(DATA_CLEAN_f(2,:)),'-r');
                hold on;
                %stem(f/1E+3,H,'-y');
                grid on;
                grid minor;
                xlim([2*(-fCent-chanHalfBW)/1E+3,2*(fCent+chanHalfBW)/1E+3]);
                xlabel('Frequency [kHz]');
                ylabel('Magnitude');
                legend('chan2');
                hold off;
            subplot(2,4,5);
                stem(f/1E+3,abs(DATA_CLEAN_f(3,:)),'-m');
                hold on;
                %stem(f/1E+3,H,'-y');
                grid on;
                grid minor;
                xlim([2*(-fCent-chanHalfBW)/1E+3,2*(fCent+chanHalfBW)/1E+3]);
                xlabel('Frequency [kHz]');
                ylabel('Magnitude');
                legend('chan3');
                hold off;
            subplot(2,4,6);
                stem(f/1E+3,abs(DATA_CLEAN_f(4,:)),'-g');
                hold on;
                %stem(f/1E+3,H,'-y');
                grid on;
                grid minor;
                xlim([2*(-fCent-chanHalfBW)/1E+3,2*(fCent+chanHalfBW)/1E+3]);
                xlabel('Frequency [kHz]');
                ylabel('Magnitude');
                legend('chan4');
                hold off;
                
            subplot(2,4,3);
                plot(t*1E+6,DATA_CLEAN_t(1,:),'-b');
                hold on;
                line([t(1)*1E+6,t(end)*1E+6],[THD,THD],'Color',[0,0,0]);
                grid on;
                grid minor;
                xlabel('Time [\mus]');
                ylabel('Amplitude');
                legend('chan1');
                hold off;
            subplot(2,4,4);
                plot(t*1E+6,DATA_CLEAN_t(2,:),'-r');
                hold on;
                line([t(1)*1E+6,t(end)*1E+6],[THD,THD],'Color',[0,0,0]);
                grid on;
                grid minor;
                xlabel('Time [\mus]');
                ylabel('Amplitude');
                legend('chan2');
                hold off;
            subplot(2,4,7);
                plot(t*1E+6,DATA_CLEAN_t(3,:),'-m');
                hold on;
                line([t(1)*1E+6,t(end)*1E+6],[THD,THD],'Color',[0,0,0]);
                grid on;
                grid minor;
                xlabel('Time [\mus]');
                ylabel('Amplitude');
                legend('chan3');
                hold off;
            subplot(2,4,8);
                plot(t*1E+6,DATA_CLEAN_t(4,:),'-g');
                hold on;
                line([t(1)*1E+6,t(end)*1E+6],[THD,THD],'Color',[0,0,0]);
                grid on;
                grid minor;
                xlabel('Time [\mus]');
                ylabel('Amplitude');
                legend('chan4');
                hold off;
    end
    
    if (OS_dwellTime == false)
        pause(10);
        OS_dwellTime = true;
    end
    
    pause(dwellTime);
    
end
