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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% CONVERTING SINUSOID FROM DECIMAL TO BINARY (2's COMP %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yBin = zeros(bits,N);

yBin = DECTOBIN(y,1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% WRITING .TXT FILE %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%
%%%%% ENTITY %%%%%
%%%%%%%%%%%%%%%%%%

fileID = fopen(fileName,'w');
fprintf(fileID,'-- Author: Joshua Simmons\n');
fprintf(fileID,'--\n');
fprintf(fileID,'-- Date: %s\n', date);
fprintf(fileID,'--\n');
fprintf(fileID,'-- Description: Sinusoidal Numerically Controlled Oscillator (NCO) via Look-Up-Table (LUT).\n');
fprintf(fileID,'-- \t\t\t\tThe input x is to be connected to a mod %d unsigned up-counter (free-running).\n', xBits);
fprintf(fileID,'-- \t\t\t\tThe output y is %d bit signed (2s Comp) with the binary point inbetween bits 1 and 2.\n', yBits);
fprintf(fileID,'\n');
fprintf(fileID,'library ieee;\n');
fprintf(fileID,'use ieee.std_logic_1164.all;\n');
fprintf(fileID,'\n');
fprintf(fileID,'ENTITY %s IS\n',entityName);
fprintf(fileID,'\tPORT(\n');
%fprintf(fileID,'\t\t--clk_1\t: IN\tstd_logic; -- FOR XSG\n');
%fprintf(fileID,'\t\t--ce_1\t: IN\tstd_logic; -- FOR XSG\n');
fprintf(fileID,'\t\tx\t: IN\tstd_logic_vector(%d DOWNTO 0);\n', xBits-1);
fprintf(fileID,'\t\ty\t: OUT\tstd_logic_vector(%d DOWNTO 0)\n', yBits-1);
fprintf(fileID,'\t);\n');
fprintf(fileID,'END %s;\n', entityName);
fprintf(fileID,'\n');

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% ARCHITECTURE %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(fileID,'ARCHITECTURE %s OF %s IS\n', archName, entityName);
fprintf(fileID,'BEGIN\n');
fprintf(fileID,'\n');
fprintf(fileID,'\tWITH x SELECT\n');

fprintf(fileID,'\t\ty <= \n');

    for row = 1:xRows;
    	fprintf(fileID,'\t\t\t"');
    	
    	%%%%%%%%%%%%%%%%%%%%%
		%%%%% WRITING Y %%%%%
		%%%%%%%%%%%%%%%%%%%%%
    	for col = 1:yBits; 
            fprintf(fileID,'%d', y_bin(row,col));
    	end

    	fprintf(fileID,'" WHEN "');

    	%%%%%%%%%%%%%%%%%%%%%
		%%%%% WRITING X %%%%%
		%%%%%%%%%%%%%%%%%%%%%
        for col = 1:xBits;
            fprintf(fileID,'%d', x_bin(row,col));
        end
        
        fprintf(fileID,'",\n');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% WRITING "WHEN OTHERS" %%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(fileID,'\n\t\t\t"');
    for row = 1:yBits; 
        fprintf(fileID,'0');
    end
    fprintf(fileID,'" WHEN OTHERS;\n');

fprintf(fileID,'\n');
fprintf(fileID,'END %s;\n', archName);

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% CLOSING FILE %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
fclose(fileID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% END OF PROGRAM %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
