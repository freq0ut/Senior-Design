%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:  Joshua Simmons
% Started: September, 2015
% Status:  COMPLETED
%
% Description: Using a limited number of samples, this algorithm
%              synchronizes to the pinger.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;
clc;

addpath('C:\Users\Joshua Simmons\Desktop\Senior_Design\Senior-Design\MATLAB\Support_Functions');
%addpath('/Users/betio32/Desktop/Senior-Design/MATLAB/Support_Functions');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETER INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Global Simulation Parameters
visualization = true;

% Pinger Properties
fPing = 30E+3;   % Source freq [Hz]
tPing = 1/fPing; % Source period [s]
tOn = 1E-3;      % Pulse On Duration [s]
PRT = 10E-3;     % Pulse-Repetitive-Period [s]

% ADC
fADC = 1800E+3;  % Sample freq [Hz]
tADC = 1/fADC;   % Sample period [s]
frameSize = 512; % Samples per frame

% Processor Properties
headCount = 0;
PRT_Array = 0.8*PRT * ones(1,10);
tError = 0;
THD = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTRUCTING INPUT SIGNAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = 0:tADC:tOn; % Time array [s]

y = sin(2*pi*fPing*t); % One pulse without zeros
y = padarray(y,[0,round((PRT-tOn)/tADC)],0,'post'); % One pulse with zeros added
y = repmat(y,[1,20]); % Copying multiple pulses
y = padarray(y,[0,round(rand()*10*frameSize)],0,'pre'); % Random start

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN SIGNAL PROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iGlobal = 1;
tGlobalLast = -1;
while ( iGlobal+frameSize-1 <= length(y) )
    
    yPartial = y(iGlobal:iGlobal+frameSize-1);
    tHeads = BREAK_THRESHOLD(yPartial,THD,'LR') * tADC;
    
    if ( tHeads > tPing/4 && tHeads < (tADC*frameSize) - tPing/4 )   
        
        % Counting heads for median array
        headCount = headCount + 1;
        
        % The start of the present frame globally
        tGlobal = tADC*(iGlobal-1);
        
        % Adjusting the PRT
        if (headCount > 1 && headCount < 12)
           PRT_Array(mod(headCount,10)+1) = tGlobal - tGlobalLast;
           PRT = median(PRT_Array);
        end
                
        if (headCount < 12) % Delaying by PRT
            tGlobalLast = tGlobal;
            tGlobal = tGlobal + PRT;
        else % Delaying so that the heads are in the center of the frame
            tError = tHeads - tADC*(frameSize/2);
            tGlobal = tGlobal + PRT + tError + tPing/4;   
        end
        
        if ( visualization == true && headCount == 13)
            figure(1);
                t = 0:tADC:(frameSize-1)*tADC; % Reconstructing time array
                yPartial = y(iGlobal:iGlobal+frameSize-1);
                stem(t*1E+6,yPartial,'-r')
                grid on;
                grid minor;
                xlabel('Time [\mus]');
                ylabel('Voltage [V]');
                title('y_{Partial}');
        end
        
        % Where we will be next globally
        iGlobal = round(tGlobal/tADC) - frameSize;
        
    end   
    
    iGlobal = iGlobal+frameSize;
end
