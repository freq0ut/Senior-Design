%   Joshua Simmons
%
%   January 13, 2016
%
%   Analysis of Microchip FFT Twiddle Factors provided in their sample code.

close all;
clear all;
clc;

addpath('/Users/betio32/Documents/myGitHub/Senior-Design/MATLAB/Support_Functions/FFT/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% READ IN DATA FROM .TXT FILE %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileID = fopen('Twiddle_Factors_256.txt','r');% Open file.

    Data_File = textscan(fileID,'%s %s');% Read in .TXT file into memory as a cell array.
    
    N_TFS = size(Data_File{1,1},1);% Total number of Twiddle Factors.
    Twids_Microchip = zeros(1,N_TFS);% Declare and initialize Microchip Twiddle Factor array.

    for i = 0:N_TFS-1;
        % In-Phase
        tempCell = textscan(Data_File{1,1}{i+1}(3:6),'%c');% Extract a specific cell in the cell array.
        tempStringHex = char(tempCell);% Convert the cell to a string of characters (hexadecimal)
        tempBin = HEX_2_BIN(tempStringHex);% Convert the string of characters (hexadecimal) to integer (binary)
        Twids_Microchip(i+1) = BIN_2_DEC(tempBin);% Convert the integer binary (2s Comp) to integer decimal 

        % Quadrature
        tempCell = textscan(Data_File{1,2}{i+1}(3:6),'%c');
        tempStringHex = char(tempCell);
        tempBin = HEX_2_BIN(tempStringHex);
        Twids_Microchip(i+1) = Twids_Microchip(i+1) + 1i*BIN_2_DEC(tempBin);
    end

    Twids_Microchip = Twids_Microchip / 2^15;
    
fclose(fileID);% Close file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% MATLAB TWIDDLE FACTORS FOR COMPARISON %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

deltaTheta = 2*pi/(2*N_TFS);

Twids_MATLAB = zeros(1,N_TFS);

for i=0:N_TFS-1;
    Twids_MATLAB(i+1) = exp(-1i*i*deltaTheta);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% PLOTTING %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig1 = figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(2,1,1);
        hold on;
        plot(real(Twids_Microchip),imag(Twids_Microchip),'-b');
        plot(real(Twids_MATLAB),imag(Twids_MATLAB),'-r');
        grid on;
        grid minor;
        xlabel('Real');
        ylabel('Imaginary');
        legend({'MICROCHIP','MATLAB'});
        titleString = sprintf('Twiddle Factors');
        title(titleString);
        hold off;
    subplot(2,1,2);
        hold on;
        plot(real(Twids_Microchip)-real(Twids_MATLAB),'-b');
        plot(imag(Twids_Microchip)-imag(Twids_MATLAB),'-r');
        grid on;
        grid minor;
        xlabel('Index');
        ylabel('Error');
        legend({'Real','Imag'});
        titleString = sprintf('Twiddle Factor Error');
        title(titleString);
        hold off;
