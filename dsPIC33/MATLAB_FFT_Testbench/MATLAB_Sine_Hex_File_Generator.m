% Author: Joshua Simmons
%
% Date: January 15, 2016
%
% Description: Creates a simulated sampled input file for use with the Microchip FFT example code.

close all;
clear all;
clc;

addpath('/Users/betio32/Documents/GitHub/myPrograms/myMATLAB/my_MATLAB_Functions');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% USER INPUT %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bits = 16;% Signed (2s Comp) with binary point between bits 1 and 2
fS = 1e6;% Sample Frequency [Hz]
N = 256;% Total number of samples

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% GENERATING SINUSOID %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ts = 1/fS;

t = Ts*(0:N-1);

y1 = cos(2*pi*40e3*t);
y2 = cos(2*pi*40e3*(t-5e-6));

for i=1:N;
    if (y1(i) == 1)
        y1(i) = 1-2^-15;
    end
    
    if (y2(i) == 1)
        y2(i) = 1-2^-15;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% CONVERTING SINUSOID FROM DECIMAL TO BINARY (2's COMP) TO HEX %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y1Bin = DECTOBIN(y1,16,1);
y1Hex = BINTOHEX(y1Bin);

y2Bin = DECTOBIN(y2,16,1);
y2Hex = BINTOHEX(y2Bin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% WRITING .TXT FILE %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

fileID = fopen('Chan2_Sampled.txt','w');

hexs = size(y1Hex,1);

for col=1:2*N;

    if (col <= N)
        fprintf(fileID,'0x');
        
        for row=1:hexs;
            fprintf(fileID,'%c',char(y2Hex(row,col)));
        end
        
        fprintf(fileID,', ');
    else
        fprintf(fileID,'0x0000, ');
    end
    
    if (mod(col,5) == 0)
        fprintf(fileID,'\n');
    end
    
end

fclose(fileID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% END OF PROGRAM %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
