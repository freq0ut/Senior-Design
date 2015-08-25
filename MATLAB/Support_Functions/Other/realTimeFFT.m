close all;
clear all;
clc;

% Simulation Parameters
tStop = 60;

% FFT Parameters
fS = 48e3;              % Sampling Frequency [Hz]
tS = 1/fS;              % Sampling Time [s]
N0 = 1024;              % Total number of data points (power of 2)
t0 = N0*tS;             % Truncation Time Interval [s]
f0 = 1/t0;              % Frequency Resolution [Hz]
t  = 0:tS:tS*(N0-1);    % Time Array
f  = f0*(-N0/2:N0/2-1); % Frequency Array (double-sided FFT)

% Ideal Filter Design
fCent = 0;  % Center Frequenc y [Hz]
fChan = fS; % Channel Half Width [Hz]
H = zeros(N0,1);
for k = 1:N0
    if ( abs(f(k)) > fCent-fChan && abs(f(k)) < fCent+fChan)
        H(k) = 1;
    end
end

% Type "audiodevinfo" for more information
chan1 = audiorecorder(fS,16,1,1); % Sampling Rate, #Bits, #Channels, DeviceID
chan2 = audiorecorder(fS,16,1,2); % Sampling Rate, #Bits, #Channels, DeviceID
chan3 = audiorecorder(fS,16,1,3); % Sampling Rate, #Bits, #Channels, DeviceID
chan4 = audiorecorder(fS,16,1,4); % Sampling Rate, #Bits, #Channels, DeviceID

% Main Recording and Visualization Loop
frameTotal=round(tStop/(N0*tS));
for frameCounter=1:frameTotal
    % Recording Audio
    recordblocking(chan1, N0*tS); % This isn't wanting to work
    recordblocking(chan2, N0*tS);
    recordblocking(chan3, N0*tS);
    recordblocking(chan4, N0*tS); % This isn't wanting to work
    
    % Device Error: Unanticipated host error???
    
    % Retrieving Recorded Audio Data    
    chan1Data = getaudiodata(chan1);
    chan2Data = getaudiodata(chan2);
    chan3Data = getaudiodata(chan3);
    chan4Data = getaudiodata(chan4);

    % Calculating FFT (double-sided)
    CHAN1DATA = fftshift(fft(chan1Data))/N0;
    CHAN2DATA = fftshift(fft(chan2Data))/N0;
    CHAN3DATA = fftshift(fft(chan3Data))/N0;
    CHAN4DATA = fftshift(fft(chan4Data))/N0;
    
    % Filtering
    CHAN1DATA_F = H.*CHAN1DATA;
    CHAN2DATA_F = H.*CHAN2DATA;
    CHAN3DATA_F = H.*CHAN3DATA;
    CHAN4DATA_F = H.*CHAN4DATA;    
    
    % Transforming back into time-domain
    chan1Data_f = ifftshift(ifft(CHAN1DATA_F))*N0;
    chan2Data_f = ifftshift(ifft(CHAN2DATA_F))*N0;
    chan3Data_f = ifftshift(ifft(CHAN3DATA_F))*N0;
    chan4Data_f = ifftshift(ifft(CHAN4DATA_F))*N0;

    % Visualization
    figure(1);
        subplot(2,2,1);
            hold on;
            plot(t*1000,chan1Data_f,'b');
            plot(t*1000,chan2Data_f,'r');
            plot(t*1000,chan3Data_f,'m');
            plot(t*1000,chan4Data_f,'g');
            hold off;
            grid;
            grid minor;
            xlabel('Time [ms]');
            xlim([0, N0*tS*1000]);
            ylabel('Amplitude [?]');
            ylim([-0.5, 0.5]);
            title({'Filtered Audio Signal','TIME DOMAIN'});
        subplot(2,2,2);
            hold on;
            xcorr(chan1Data_f,chan2Data_f,'coeff');
            xcorr(chan1Data_f,chan3Data_f,'coeff');
            xcorr(chan1Data_f,chan4Data_f,'coeff');
            legend({'Chan12','Chan13','Chan14'});
            title('Cross-Correlation');
            hold off;
        subplot(2,2,3);
            hold on;
            plot(f/1000,20*log10(abs(CHAN1DATA_F)),'b');
            plot(f/1000,20*log10(abs(CHAN2DATA_F)),'r');
            plot(f/1000,20*log10(abs(CHAN3DATA_F)),'m');
            plot(f/1000,20*log10(abs(CHAN4DATA_F)),'g');
            hold off;
            grid;
            grid minor;
            xlabel('Frequency [kHz]');
            xlim([-fS/2000, fS/2000]);
            ylabel('Magnitude [dB]');
            ylim([-140, 0]);
            legend({'Chan1','Chan2','Chan3','Chan4'});
            title({'Filtered Audio Signal','FREQUENCY DOMAIN'});
        subplot(2,2,4);
            compass([0,1],[0,0]);
            view([90 -90]);
            title({'Estimated Direction of Target Signal'});            
end