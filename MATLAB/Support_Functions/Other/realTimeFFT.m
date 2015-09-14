close all;
clear all;
clc;

% Simulation Parameters
tStop = 10*60;
oneShotWindowSizing = false;

% FFT Parameters
fS = 192E+3;            % Sampling Frequency [Hz]
tS = 1/fS;              % Sampling Time [s]
N0 = 4096;              % Total number of data points (power of 2)
t0 = N0*tS;             % Truncation Time Interval [s]
f0 = 1/t0;              % Frequency Resolution [Hz]
t  = 0:tS:tS*(N0-1);    % Time Array
f  = f0*(-N0/2:N0/2-1); % Frequency Array (double-sided FFT)
f  = f';

% Ideal Bandpass Filter
fCent = 0;  % Center Frequency [Hz]
fChan = fS; % Channel Half Width [Hz]
H = zeros(N0,1);
for i = 1:N0
    if ( abs(f(i)) >= fCent-fChan && abs(f(i)) <= fCent+fChan)
        H(i) = 1;
    end
end

% Type "audiodevinfo" for more information
chan13 = audiorecorder(fS,24,1,2); % Sampling Rate, #Bits, #Channels, DeviceID
% chan24 = audiorecorder(fS,24,2,2); % Sampling Rate, #Bits, #Channels, DeviceID

% Main Recording and Visualization Loop
frameTotal = ceil(tStop/(N0*tS));

for frameCounter = 1:frameTotal;
    
    % Recording Audio (synchronously)
    recordblocking(chan13, N0*tS);
%     recordblocking(chan24, N0*tS);
    
%     % Recording Audio (asynchronously)
%     record(chan13);
%     record(chan24);
%     pause(1E-3);
%     
%     stop(chan13);
%     stop(chan24);
    
    % Retrieving Recorded Audio Data    
    DATA_RAW_t_1 = getaudiodata(chan13);
%     DATA_RAW_t_3 = DATA_RAW_t_1(:,2);
%     DATA_RAW_t_1(:,2) = [];
%     
%     DATA_RAW_t_2 = getaudiodata(chan24);
%     DATA_RAW_t_4 = DATA_RAW_t_2(:,2);
%     DATA_RAW_t_2(:,2) = [];
    
    % Truncating data to length N0
    DATA_RAW_t_1(N0+1:end) = [];
%     DATA_RAW_t_2(N0+1:end) = [];
%     DATA_RAW_t_3(N0+1:end) = [];
%     DATA_RAW_t_4(N0+1:end) = [];

    % Calculating FFT (double-sided)
    DATA_RAW_f_1 = fftshift(fft(DATA_RAW_t_1)) / N0;
%     DATA_RAW_f_2 = fftshift(fft(DATA_RAW_t_2)) / N0;
%     DATA_RAW_f_3 = fftshift(fft(DATA_RAW_t_3)) / N0;
%     DATA_RAW_f_4 = fftshift(fft(DATA_RAW_t_4)) / N0;
%     
    % Filtering
%     DATA_CLEAN_f_1 = H.*DATA_RAW_f_1;
%     DATA_CLEAN_f_3 = H.*DATA_RAW_f_3;
    
    % Transforming back into time-domain
%     DATA_CLEAN_t_1 = ifft(ifftshift(DATA_CLEAN_f_1)) * N0;
%     DATA_CLEAN_t_3 = ifft(ifftshift(DATA_CLEAN_f_3)) * N0;

    % Calculating signal power
    power1 = DATA_RAW_f_1 .* conj(DATA_RAW_f_1);
%     power2 = DATA_RAW_f_2 .* conj(DATA_RAW_f_2);
%     power3 = DATA_RAW_f_3 .* conj(DATA_RAW_f_3);
%     power4 = DATA_RAW_f_4 .* conj(DATA_RAW_f_4);

    % Visualization
    figure(1);
%         subplot(2,2,1);
            plot(f/1E+3,10*log10(power1*1E+3),'b');
            grid;
            grid minor;
            xlabel('Frequency [kHz]');
            xlim([-20, 20]);
            ylabel('Power [dBm]');
            ylim([-60, 0]);
            legend({'Chan1'});
%         subplot(2,2,2);
%             plot(f/1E+3,10*log10(power2*1E+3),'r');
%             grid;
%             grid minor;
%             xlabel('Frequency [kHz]');
%             xlim([-20, 20]);
%             ylabel('Power [dBm]');
%             ylim([-60, 0]);
%             legend({'Chan2'});
%         subplot(2,2,3);
%             plot(f/1E+3,10*log10(power3*1E+3),'m');
%             grid;
%             grid minor;
%             xlabel('Frequency [kHz]');
%             xlim([-20, 20]);
%             ylabel('Power [dBm]');
%             ylim([-60, 0]);
%             legend({'Chan3'});
%         subplot(2,2,4);
%             plot(f/1E+3,10*log10(power4*1E+3),'g');
%             grid;
%             grid minor;
%             xlabel('Frequency [kHz]');
%             xlim([-20, 20]);
%             ylabel('Power [dBm]');
%             ylim([-60, 0]);
%             legend({'Chan4'});
            
        if (oneShotWindowSizing == true)
%         else
            oneShotWindowSizing = false;
%             pause(5);
        end
end