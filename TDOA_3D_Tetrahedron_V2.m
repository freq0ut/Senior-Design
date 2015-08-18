%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Joshua Simmons
% August 17, 2015
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
%   S  = (xS,yS,zS)
%
% Sequence of Events
% 1. Initialization of parameters.
% 2. Random number generator computes time delays for chan2, chan3 and
%    chan4.
% 3. Input signals are constructed with the aforementioned time delays
%    and white Gaussian noise is added.
% 4. Cross-correlations (XC) are computed for chan2, chan3, and chan4
%    using chan1 as the reference.
% 5. The maximum y-coordinate of each XC is found and the corresponding
%    x-coordinate is multiplied by the sample time. This is the
%    estimated time delay.
% 6. The time delays are plugged into formulas to find the grid
%    coordinates of the source.
% 7. Due to the random number generator, sometimes impossible source
%    locations arise. These results turn out to be complex with an
%    imaginary component. Because of this a filter is applied so that
%    only real results pass thru. If the source location is real then
%    the trialCounter increments else it does not.
% 8. Once the grid coordinates of the source are found, the horizontal and 
%    vertical azimuths are computed.
% 9. Results are visualized.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;
clc;

% Global Simulation Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trialTotal = 1;
fig1_On = true;  % Time domain signals and cross-correlation functions
fig2_On = true;  % Compass plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Source Properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SNR  = 10;        % Signal to Noise Ratio [dB]
fSce = 30e3;      % Source freq [Hz]
Tsce = 1/fSce;    % Source period [s]
vP   = 1482;      % Propagation Velocity [m/s]
lambda = vP/fSce; % Wavelength [m]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Hydrophone Properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = (1/3)*lambda; % Hydrophone spacing [m]
d = D / sqrt(2);  % For coordinates of hydrophones [m]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Sampler Properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fS = 1e6;  % Sample freq [Hz]
tS = 1/fS; % Sample period [s]

% POWER OF 2 STRONGLY PREFERRED TO TAKE ADVANTAGE OF FFT ACCELERATIONS
N0 = 2^7; % #Samples per frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% START MAIN LOOP
for trialCount = 1:trialTotal;    
    % Constructing simulated input signals
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Using a RNG to randomize the actual source location
    for iter = 1:3;
        rNum = rand() * 1e3; % Between 0 and 1000
    
        if ( rNum > 1 )
            S_Act(iter) = rNum; % [xS_Act, yS_Act, zS_Act]
        end
    end

    % Determining the actual time delays
    R1_Act = sqrt( (S_Act(1)-d)^2 + (S_Act(2)  )^2 + (S_Act(3)  )^2 );
    R2_Act = sqrt( (S_Act(1)  )^2 + (S_Act(2)-d)^2 + (S_Act(3)  )^2 );
    R3_Act = sqrt( (S_Act(1)  )^2 + (S_Act(2)  )^2 + (S_Act(3)-d)^2 );
    R4_Act = sqrt( (S_Act(1)-d)^2 + (S_Act(2)-d)^2 + (S_Act(3)-d)^2 );
    
    td2_Act = (R2_Act-R1_Act) / vP;
    td3_Act = (R3_Act-R1_Act) / vP;
    td4_Act = (R4_Act-R1_Act) / vP;
    
    % Constructing the input signals from the time delays
    t = 0:tS:(N0-1)*tS;       % Time Array [s]
    chan1 = cos(2*pi*fSce*t); % Reference Signal
    chan2 = cos(2*pi*fSce*(t+td2_Act));
    chan3 = cos(2*pi*fSce*(t+td3_Act));
    chan4 = cos(2*pi*fSce*(t+td4_Act));

    % Adding Gaussian Noise for increased realism
    chan1 = awgn(chan1,SNR);
    chan2 = awgn(chan2,SNR);
    chan3 = awgn(chan3,SNR);
    chan4 = awgn(chan4,SNR);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % Determining the estimated time delays from cross-correlation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [XC12, XC12lags] = xcorr(chan1,chan2,'coeff');
    [XC12yMax, XC12yMaxK] = max(XC12);
    td2_Est  = XC12lags(XC12yMaxK)*tS;

    [XC13, XC13lags] = xcorr(chan1,chan3,'coeff');
    [XC13yMax, XC13yMaxK] = max(XC13);
    td3_Est  = XC13lags(XC13yMaxK)*tS;

    [XC14, XC14lags] = xcorr(chan1,chan4,'coeff');
    [XC14yMax, XC14yMaxK] = max(XC14);
    td4_Est  = XC14lags(XC14yMaxK)*tS;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % Determing the location and azimuth to the source
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Calculating estimated source location
    r2_Est = vP*td2_Est;
    r3_Est = vP*td3_Est;
    r4_Est = vP*td4_Est;
    
    R1_Est = ( (r2_Est)^2+(r3_Est)^2+(r4_Est)^2+(4*d)^2 ) / ...
        ( -2*(r2_Est+r3_Est+r4_Est) );
    
    R2_Est = R1_Est + r2_Est;
    R3_Est = R1_Est + r3_Est;
    R4_Est = R1_Est + r4_Est;
    
    xS_Est = (-R1_Est^2 + R2_Est^2 + R3_Est^2 - R4_Est^2 + 2*d^2) / (4*d);
    yS_Est = ( R1_Est^2 - R2_Est^2 + R3_Est^2 - R4_Est^2 + 2*d^2) / (4*d);
    zS_Est = ( R1_Est^2 + R2_Est^2 - R3_Est^2 - R4_Est^2 + 2*d^2) / (4*d);

    % Converting azimuths from radians to degrees and wrapping them so
    % that they fall between 0 and 360 degrees
    azimuthH_Act = wrapTo2Pi(atan2(S_Act(2),S_Act(1))) * (180/pi);
    azimuthV_Act = wrapTo2Pi(atan2(S_Act(3),S_Act(1))) * (180/pi);
        
    azimuthH_Est = wrapTo2Pi(atan2(yS_Est,xS_Est)) * (180/pi);
    azimuthV_Est = wrapTo2Pi(atan2(zS_Est,xS_Est)) * (180/pi);

    % Visualization
        
    stringTrials = sprintf('Trial %0.0f / %0.0f', trialCount, trialTotal);
        
    if (fig1_On == true)
    
        % Stemming time domain (superimposed) and cross-correlations
        figure(1)
            subplot(2,2,1);
                hold on;
                stem(t*1e6,chan1,'-b');
                stem(t*1e6,chan2,'-r');
                stem(t*1e6,chan3,'-m');
                stem(t*1e6,chan4,'-g');
                line([Tsce/2*1e6, Tsce/2*1e6],[-1.5, 1.5]);
                hold off;
                xlabel('Time [\mus]');
                ylabel('Amplitude');
                legend({'Chan1','Chan2','Chan3','Chan4'});
                string111 = sprintf('f_{samp} = %0.0f [MHz]', fS/1e6);
                string112 = sprintf('SNR = %0.0f [dB]', SNR);
                title({string111,string112});
            subplot(2,2,2);
                hold on;
                stem(XC12lags*tS*1e6,XC12,'r');
                plot(td2_Act*1e6,XC12yMax,'b.','MarkerSize',20);
                plot(td2_Est*1e6,XC12yMax,'k.','MarkerSize',20);
                hold off;
                string121 = sprintf('td2_{Act} = %f [us]', td2_Act*1e6);
                string122 = sprintf('td2_{Est} = %f [us]', td2_Est*1e6);
                legend({'',string121,string122});
                title('XC_{12}');
                xlabel('Time [\mus]');
            subplot(2,2,3);
                hold on;
                stem(XC13lags*tS*1e6,XC13,'m');
                plot(td3_Act*1e6,XC13yMax,'b.','MarkerSize',20);
                plot(td3_Est*1e6,XC13yMax,'k.','MarkerSize',20);
                hold off;
                string131 = sprintf('td3_{Act} = %f [us]', td3_Act*1e6);
                string132 = sprintf('td3_{Est} = %f [us]', td3_Est*1e6);
                legend({'',string131,string132});
                title('XC_{13}');
                xlabel('Time [\mus]');
            subplot(2,2,4);
                hold on;
                stem(XC14lags*tS*1e6,XC14,'g');
                plot(td4_Act*1e6,XC14yMax,'b.','MarkerSize',20);
                plot(td4_Est*1e6,XC14yMax,'k.','MarkerSize',20);
                hold off;
                string141 = sprintf('td4_{Act} = %f [us]', td4_Act*1e6);
                string142 = sprintf('td4_{Est} = %f [us]', td4_Est*1e6);
                legend({'',string141,string142});
                title('XC_{14}');
                xlabel('Time [\mus]');
    end
        
    if (fig2_On == true)
            
        % Comparint actual and estimated source azimuths
        figure(2)
            subplot(2,2,1);
                compass([0 S_Act(1)],[0 S_Act(2)],'-r');
                view([90, -90]);
                string211 = sprintf('Actual: %0.1f (deg)', azimuthH_Act);
                title({stringTrials,'Horizontal Azimuth',string211,''});
            subplot(2,2,3);
                compass([0 xS_Est],[0 yS_Est],'-b');
                view([90, -90]);
                string241 = sprintf('Estimated: %0.1f (deg)', azimuthH_Est);
                title({stringTrials,'Horizontal Azimuth',string241,''});
                
            subplot(2,2,2);
                compass([0 S_Act(1)],[0 S_Act(3)],'-r');
                view([90, -90]);
                string221 = sprintf('Actual: %0.1f (deg)', azimuthV_Act);
                title({stringTrials,'Vertical Azimuth',string221,''});
            subplot(2,2,4);
                compass([0 xS_Est],[0 zS_Est],'-b');
                view([90, -90]);
                string251 = sprintf('Estimated: %0.1f (deg)', azimuthV_Est);
                title({stringTrials,'Vertical Azimuth',string251,''});
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pause(2.00);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
