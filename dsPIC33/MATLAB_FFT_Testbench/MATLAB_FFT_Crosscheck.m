% Joshua Simmons
%
% January 11, 2016
%
% Crosschecking Microchip DSPic FFT Output with a square wave input.

close all;
clear all;
clc;

fileID = fopen('sigCmpx.csv','r');

    CSV_File = textscan(fileID,'%s %s');
    
    N = size(CSV_File{1,1},1)/2;
    IQ_Data = zeros(1,N);
    
    for i = 0:N-1;
        tempCell = textscan(CSV_File{1,2}{2*i+1}(3:6),'%c');
        tempStringHex = char(tempCell);
        tempBin = HEX_2_BIN(tempStringHex);
        IQ_Data(i+1) = BIN_2_DEC(tempBin);
        
        tempCell = textscan(CSV_File{1,2}{2*i+2}(3:6),'%c');
        tempStringHex = char(tempCell);
        tempBin = HEX_2_BIN(tempStringHex);
        IQ_Data(i+1) = IQ_Data(i+1) + 1i*BIN_2_DEC(tempBin);
    end

fclose(fileID);

fS = 1e3;
Ts = 1/fS;

f0 = 1/((N-1)*Ts);

f = f0*(0:N-1);
t = Ts*(0:N-1);

y = ones(1,N);
for i=0:N-1;
    if (mod(i,10) >= 5)
        y(i+1) = -1;
    end
end

Y = fft(y) / N;

fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(2,2,1);
    	stem(t,y,'-b');
    	grid on;
    	grid minor;
        xlim([t(1),t(end)]);
    	xlabel('Time [s]');
    	ylabel('Amplitude [V]');
        titleString = sprintf('TEMPORAL\nMATLAB Simulated Input Signal');
        title(titleString);
    subplot(2,2,3);
    	stem(f,abs(Y),'-r');
    	grid on;
    	grid minor;
    	xlabel('Frequency [Hz]');
    	ylabel('Magnitude [Volts]');
        titleString = sprintf('FREQUENCY\nMATLAB Simulated Input Signal');
        title(titleString);
    subplot(2,2,4);
    	stem(f,abs(IQ_Data)/2^14,'-r');
    	grid on;
    	grid minor;
    	xlabel('Frequency [Hz]');
    	ylabel('Magnitude [Volts]');
        titleString = sprintf('FREQUENCY\ndsPIC Simulated Output Data');
        title(titleString);