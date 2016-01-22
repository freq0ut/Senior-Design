% Author: Joshua Simmons
%
% Date: January 15, 2016
%
% Description: Creates a simulated sampled input file for use with the Microchip FFT example code.

close all;
clear all;
clc;

addpath('/Users/betio32/Documents/myGitHub/myPrograms/myMATLAB/my_MATLAB_Functions');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% USER INPUT %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileName = sprintf('Chan1_Sampled.txt');

bits = 16;% Signed (2s Comp) with binary point between bits 1 and 2
fS = 1e6;% Sample Frequency [Hz]
N = 256;% Total number of samples

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% GENERATING SINUSOID %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tS = 1/fS;

t = tS*(0:N-1);

y = cos(2*pi*40e3*t);

for i=1:N;
    if (y(i) == 1)
        y(i) = 1-2^-15;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% CONVERTING SINUSOID FROM DECIMAL TO BINARY (2's COMP) TO HEX %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yBin = DECTOBIN(y,16,1);
yHex = BINTOHEX(yBin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% WRITING .TXT FILE %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

fileID = fopen(fileName,'w');

hexs = size(yHex,1);

for col=1:N;
    
    fprintf(fileID,'0x');
    
    for row=1:hexs;
        fprintf(fileID,'%c',char(yHex(row,col)));
    end

    fprintf(fileID,',\n');

end

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% CLOSING FILE %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
fclose(fileID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% END OF PROGRAM %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
