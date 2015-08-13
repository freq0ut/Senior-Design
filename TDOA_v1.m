%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Joshua Simmons                                                     %
% July 2015                                                          %
% Uses Time-Difference-of-Arrival (TDOA) to determine the azimuth to %
% a 30 kHz SINE wave underwater.                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;
clc;

% Global Simulation Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trialCount = 0;
trialTotal = 20;
fig1_On = false;  % Time domain signals and cross-correlation functions
fig2_On = true;  % Compass plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Source Properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SNR  = 40;        % Signal to Noise Ratio [dB]
fSce = 30e3;      % Source freq [Hz]
Tsce = 1/fSce;    % Source period [s]
vP   = 1482;      % Propagation Velocity [m/s]
lambda = vP/fSce; % Wavelength [m]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Hydrophone Properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d = (1/3)*lambda; % Hydrophone spacing [m]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Sampler Properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fS = 1e6;  % Sample freq [Hz]
tS = 1/fS; % Sample period [s]

% POWER OF 2 STRONGLY PREFERRED TO TAKE ADVANTAGE OF FFT ACCELERATIONS
N0 = 2^10; % #Samples per frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% START MAIN LOOP
while (trialCount < trialTotal)
    
    % Constructing simulated input signals
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Using a random number generator to randomize the actual time delays
    for iter=2:4;
        rNum = rand();
    
        if (rNum <= 0.5)
            rNum = -1*rNum;
        end
    
        if (iter==2)
            td2_Act = rNum/2*Tsce;   % Actual time delay of chan2
        elseif (iter==3)
            td3_Act = rNum/2*Tsce;   % Actual time delay of chan3
        elseif (iter==4)
            td4_Act = rNum/2*Tsce;   % Actual time delay of chan4
        end
    end

    % Constructing the input signals from the randomized time delays
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
    tRef_Est = (td2_Est^2+td4_Est^2-td3_Est^2)/(2*(td3_Est-td2_Est-td4_Est));
    R1_Est   = vP*tRef_Est;
    R2_Est   = R1_Est + vP*td2_Est;
    yS_Est   = (R2_Est^2-R1_Est^2-d^2)/(2*d);
    xS_Est1  =  sqrt(R1_Est^2-yS_Est^2);
    xS_Est2  = -sqrt(R1_Est^2-yS_Est^2);
    
    % Calculating actual source location
    tRef_Act = (td2_Act^2+td4_Act^2-td3_Act^2)/(2*(td3_Act-td2_Act-td4_Act));
    R1_Act   = vP*tRef_Act;
    R2_Act   = R1_Act + vP*td2_Act;
    yS_Act   = (R2_Act^2-R1_Act^2-d^2)/(2*d);
    xS_Act1  =  sqrt(R1_Act^2-yS_Act^2);
    xS_Act2  = -sqrt(R1_Act^2-yS_Act^2);

    % Checking for impossibilities
    if ( (angle(xS_Est1) == 0 || angle(xS_Est1) == pi) && (angle(xS_Act1) == 0 || angle(xS_Act1) == pi) )
        trialCount = trialCount + 1;
        
        % Converting azimuths from radians to degrees
        azimuth1_Est = atan2(yS_Est,xS_Est1) * (180/pi);
        azimuth2_Est = atan2(yS_Est,xS_Est2) * (180/pi);
    
        azimuth1_Act = atan2(yS_Act,xS_Act1) * (180/pi);
        azimuth2_Act = atan2(yS_Act,xS_Act2) * (180/pi);
        
        
        % Handling negative azimuths
        if (azimuth1_Est < 0)
            azimuth1_Est = azimuth1_Est + 360;
        end

        if (azimuth2_Est < 0)
            azimuth2_Est = azimuth2_Est + 360;
        end
        
        if (azimuth1_Act < 0)
            azimuth1_Act = azimuth1_Act + 360;
        end
        
        if (azimuth2_Act < 0)
            azimuth2_Act = azimuth2_Act + 360;
        end

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
                    hold off;
                    xlim([0, 4*lambda/vP*1e6]);
                    xlabel('Time [\mus]');
                    ylabel('Amplitude');
                    legend({'Chan1','Chan2','Chan3','Chan4'});
                    string111 = sprintf('f_{source} = %0.0f [kHz]          f_{samp} = %0.0f [MHz]', fSce/1e3, fS/1e6);
                    string112 = sprintf('SNR = %0.0f [dB]',         SNR);
                    string113 = sprintf('         d = %0.4f [cm]          d = %0.4f [lambda]', d*100, d/lambda);
                    string114 = sprintf('T/4 = %0.4f [us]',         Tsce/4*1e6);
                    title({string111,string112,string113,string114});
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
            
            % Estimated source azimuths
            figure(2)
                subplot(2,2,1);            
                    compass([0 xS_Act1],[0 yS_Act],'-r');
                    view([90, -90]);
                    string211 = sprintf('%0.1f (deg)', azimuth1_Act);
                    title({stringTrials,'1^{st} Actual Azimuth',string211,''});
                subplot(2,2,2);            
                    compass([0 xS_Act2],[0 yS_Act],'-r');
                    view([90, -90]);
                    string212 = sprintf('%0.1f (deg)', azimuth2_Act);
                    title({stringTrials,'2^{nd} Actual Azimuth',string212,''});
                subplot(2,2,3);            
                    compass([0 xS_Est1],[0 yS_Est],'-b');
                    view([90, -90]);
                    string221 = sprintf('%0.1f (deg)', azimuth1_Est);
                    title({stringTrials,'1^{st} Estimated Azimuth',string221,''});
                subplot(2,2,4);            
                    compass([0 xS_Est2],[0 yS_Est],'-b');
                    view([90, -90]);
                    string231 = sprintf('%0.1f (deg)', azimuth2_Est);
                    title({stringTrials,'2^{nd} Estimated Azimuth',string231,''});
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pause(0.25);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%