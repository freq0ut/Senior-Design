%   Joshua Simmons
%
%   January 11, 2016
%
%   MATLAB Testbench for cross-checking Microchip MPLABX dsPIC30F FFT
%   Example Code for a 1kHz Square Wave Signal.

close all;
clear all;
clc;

addpath('/Users/betio32/Documents/myGitHub/Senior-Design/MATLAB/Support_Functions/FFT/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% READ IN DATA FROM THE .CSV FILE %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileID = fopen('sigCmpx.csv','r');% Open file

    Data_File = textscan(fileID,'%s %s');% Read in entire CSV file into memory as a cell array
    
    N = size(Data_File{1,1},1)/2;% Total number of samples
    Y_DSPIC = zeros(1,N);% Declare and initialize IQ Data array

    for i = 0:N-1;
        % In-Phase
        tempCell = textscan(Data_File{1,2}{2*i+1}(3:6),'%c');% Extract a specific cell in the cell array
        tempStringHex = char(tempCell);% Convert the cell to a string of characters (hexadecimal)
        tempBin = HEX_2_BIN(tempStringHex);% Convert the string of characters (hexadecimal) to integer (binary)
        Y_DSPIC(i+1) = BIN_2_DEC(tempBin);% Convert the integer binary (2s Comp) to integer decimal 
        
        % Quadrature
        tempCell = textscan(Data_File{1,2}{2*i+2}(3:6),'%c');
        tempStringHex = char(tempCell);
        tempBin = HEX_2_BIN(tempStringHex);
        Y_DSPIC(i+1) = Y_DSPIC(i+1) + 1i*BIN_2_DEC(tempBin);
    end

fclose(fileID);% Close file

% Scaling data
Y_DSPIC = Y_DSPIC / 2^14;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% TAKING THE INVERSE FFT OF THE MPLABX DSPIC30F DATA %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y_DSPIC = ifft(Y_DSPIC) * N;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% CREATING MATLAB SIMULATED SIGNAL FOR COMPARISON %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fS = 10e3;% Sampling frequency [Hz]
Ts = 1/fS;% Sampling interval [s]

f0 = fS/(N-1);% Frequency resolution [Hz]

f = f0*(0:N-1);% Frequency Array (single-sided)
t = Ts*(0:N-1);% Time Array

% Creating a 1kHz square wave (50% duty cycle) with 10 samples per period
% with ZERO DC offset.
y_MATLAB = ones(1,N);
for i=0:N-1;
    if (mod(i,10) >= 5)
        y_MATLAB(i+1) = -1;
    end
end

Y_MATLAB = fft(y_MATLAB) / N;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% PLOTTING %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(2,2,1);
    	stem(t,y_MATLAB,'-b');
    	grid on;
    	grid minor;
        xlim([t(1),t(end)]);
        ylim([-1.1,1.1]);
    	xlabel('Time [s]');
    	ylabel('Amplitude [V]');
        titleString = sprintf('TEMPORAL\nMATLAB');
        title(titleString);
    subplot(2,2,2);
    	stem(t,y_DSPIC,'-b');
    	grid on;
    	grid minor;
        xlim([t(1),t(end)]);
        ylim([-1.1,1.1]);
    	xlabel('Time [s]');
    	ylabel('Amplitude [V]');
        titleString = sprintf('TEMPORAL\ndsPIC');
        title(titleString);
    subplot(2,2,3);
    	stem(f,abs(Y_MATLAB),'-r');
    	grid on;
    	grid minor;
    	xlabel('Frequency [Hz]');
    	ylabel('Magnitude [Volts]');
        titleString = sprintf('FREQUENCY\nMATLAB');
        title(titleString);
    subplot(2,2,4);
    	stem(f,abs(Y_DSPIC),'-r');
    	grid on;
    	grid minor;
    	xlabel('Frequency [Hz]');
    	ylabel('Magnitude [Volts]');
        titleString = sprintf('FREQUENCY\ndsPIC');
        title(titleString);
