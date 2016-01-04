close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% USER INPUT %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ADC
N = 2^10;% Total number of samples
fS = 1e6;% ADC sampling frequency [Hz]

% Ultrasonic Acoustic Bacon
fPing = 20e3;% Pinger frequency [Hz]
SNR_dB = 10;% Received signal signal-to-noise-ratio [dB]
t0 = 0.4/fPing;% Time delay [s]

% Digital BPF
BPF_fC = fPing;% Digital BPF center frequency [Hz]
BPF_W = 0.2*fPing;% Digital BPF bandwidth [Hz]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% CREATING SIMULATED SIGNALS %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Creating time array
Ts = 1/fS;
t = 0:N-1;
t = Ts*t;

% Creating frequency array
T0 = (N-1)*Ts;
f0 = 1/T0;
f = -N:N-1;
f = f0/2*f;

% Transmitted waveforms
s1 = cos(2*pi*fPing*t);
s2 = cos(2*pi*fPing*(t-t0));

% Received waveforms
r1 = awgn(s1,SNR_dB);
r2 = awgn(s2,SNR_dB);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% SIGNAL PROCESSING %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Non-circular convolution
R1 = fftshift(fft(r1,2*N));
R2 = fftshift(fft(fliplr(r2),2*N));% r2 is time reversed!

% Cross-correlation (frequency domain)
XC_f = R1.*R2 / (N^2);

% Bandpass Filtering
for i=1:length(f);
    if ( abs(f(i)) < BPF_fC - BPF_W/2 || abs(f(i)) > BPF_fC + BPF_W/2)
        XC_f(i) = 0;
    end
end

% Cross-correlation (time domain)
XC_t = ifft(ifftshift(XC_f),2*N) * (2*N);

% Estimating the time delay
[~,it0] = max(XC_t);
tDelay = -Ts*(-N+it0);% Corrected (should be -N-1+it0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% PLOTTING %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(2,2,1);
        stem(t*1e3,r1,'-b');
        grid on;
        grid minor;
        xlim([t(1)*1e3,t(end)*1e3]);
        xlabel('Time [ms]');
        ylabel('Amplitude [V]');
        title('Chan_1 Time Signal');
    subplot(2,2,3);
        stem(t*1e3,r2,'-r');
        grid on;
        grid minor;
        xlim([t(1)*1e3,t(end)*1e3]);
        xlabel('Time [ms]');
        ylabel('Amplitude [V]');
        title('Chan_2 Time Signal');
    subplot(2,2,2);
        stem(f/1e3,abs(XC_f),'-g');
        grid on;
        grid minor;
        xlim([-2*fPing/1e3,2*fPing/1e3]);
        xlabel('Frequency [kHz]');
        ylabel('Magnitude [?]');
        title('Post-Filtered XC (Frequency Domain)');
    subplot(2,2,4);
        stem(Ts*(-N:N-1)*1e6,XC_t,'-M');
        grid on;
        grid minor;
        xlim([-2*t0*1e6,2*t0*1e6]);
        xlabel('Time [\mus]');
        ylabel('Correlation Coefficient (\rho)');
        title('XC (Time Domain)');
