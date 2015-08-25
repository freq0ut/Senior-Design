close all;
clear all;
clc;

% Audio Recording Info
fileName = sprintf('/Users/betio32/Desktop/AudioRecoding.wav');
tStop = 10;

AR = dsp.AudioRecorder('OutputNumOverrunSamples',true);                    
AW = dsp.AudioFileWriter(fileName,'FileFormat', 'WAV');

disp('Speak into microphone now');
tic;
while toc < tStop,
  [audioIn,nOverrun] = step(AR);
  step(AW,audioIn);
  if nOverrun > 0
    fprintf('Audio recorder queue was overrun by %d samples\n', nOverrun);
  end
end
release(AR);
release(AW);
disp('Recording complete');

g = audioread(fileName);

% FFT Parameters
fS = 44100;        % Sampling Frequency [Hz]
tS = 1/fS;         % Sampling Time [s]
t0 = length(g)*tS; % Truncation Time Interval [s]
f0 = 1/t0;         % Frequency Resolution [Hz]
N0 = length(g);    % Total number of data points

t = 0:tS:tS*(length(g)-1);
f = f0*(-N0/2:N0/2-1);

G = fftshift(fft(g)) / N0;

% Calculating the signal power in both the time-domain and the frequency-domain. 
% pg  = trapz(t,g.*conj(g)) / t0;
% pG  = trapz(f,G.*conj(G)) / f0;

% Visualization
figure(1);
    subplot(2,1,1);
        plot(t,g,'b');
        grid;
        grid minor;
        xlabel('Time [seconds]');
        xlim([t(1), t(end)]);
        ylabel('? [?]');
        title({fileName;'TIME DOMAIN'});
    subplot(2,1,2);
        plot(f,abs(G),'b');
        grid;
        grid minor;
        xlabel('Frequency [Hertz]');
        xlim([f(1), f(end)]);
        ylabel('Magnitude');
        title({fileName;'FREQUENCY DOMAIN'});
        
%sound(g, fS);