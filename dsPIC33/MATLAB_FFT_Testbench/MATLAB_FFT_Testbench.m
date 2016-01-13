% Joshua Simmons
%
% January 11, 2016
%
% MATLAB Testbench for cross-checking Microchip MPLABX dsPIC30F FFT Example Code for a Square Wave Signal.

close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% READ IN MPLABX CSV HEX FILE %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileID = fopen('sigCmpx.csv','r');% Open file

    CSV_File = textscan(fileID,'%s %s');% Read in entire CSV file into memory as a cell array
    
    N = size(CSV_File{1,1},1)/2;% Total number of samples
    IQ_Data = zeros(1,N);% Declare and initialize IQ Data array

    for i = 0:N-1;
        % In-Phase
        tempCell = textscan(CSV_File{1,2}{2*i+1}(3:6),'%c');% Extract a specific cell in the cell array
        tempStringHex = char(tempCell);% Convert the cell to a string of characters (hexadecimal)
        tempBin = HEX_2_BIN(tempStringHex);% Convert the string of characters (hexadecimal) to integer (binary)
        IQ_Data(i+1) = BIN_2_DEC(tempBin);% Convert the integer binary (2s Comp) to integer decimal 
        
        % Quadrature
        tempCell = textscan(CSV_File{1,2}{2*i+2}(3:6),'%c');
        tempStringHex = char(tempCell);
        tempBin = HEX_2_BIN(tempStringHex);
        IQ_Data(i+1) = IQ_Data(i+1) + 1i*BIN_2_DEC(tempBin);
    end

fclose(fileID);% Close file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% CREATING MATLAB SIMULATED SIGNAL FOR COMPARISON %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fS = 1e3;% Sampling frequency [Hz]
Ts = 1/fS;% Sampling interval [s]

f0 = 1/((N-1)*Ts);% Frequency resolution [Hz]

f = f0*(0:N-1);% Frequency Array (single-sided)
t = Ts*(0:N-1);% Time Array

% Creating the temporal square wave with 50% duty cycle 10 samples per period
y = ones(1,N);
for i=0:N-1;
    if (mod(i,10) >= 5)
        y(i+1) = -1;
    end
end

Y = fft(y) / N;% Taking the Fast-Fourier-Transform of the square wave

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% PLOTTING %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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