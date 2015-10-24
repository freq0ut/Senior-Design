% Author: Joshua Simmons
%
% Date: October 22, 2015
%
% Desciption:   Selects the R and C values as well as plot the freq response for
%               the Multiple Feedback Bandpass Filter.
%
%               The OP AMP that you select needs to have an open-loop gain of at
%               least 10 times above the peak amplitude of the frequency response.
%
%               In other words A0 >> H*Q

clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% USER INPUT %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f0 = 20.0E+3;% Center frequency of passband

% AT f0 THE GAIN WILL BE THE PRODUCT OF H AND Q
H = 1;% Gain at f0
Q = 10;% Quality factor

% You select C1. Keep an eye on the calculated component
% values of the resistors and capacitors so that they
% are realizeable. 
C1 = 1.11E-9;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% CALCULATING COMPONENT VALUES %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w0 = 2*pi*f0;

zeta = 1/(2*Q); % Damping ratio
alpha = w0*zeta;% Normalized damping ratio

k = C1*w0;% Just a constant to simply math expressions
C2 = C1;
R1 = 1/(H*k);
R2 = 1/((2*Q-H)*k);
R3 = 2*Q/k;

% Printing values of components
componentValueString = sprintf('f0 = %2.2f [kHz]', 0.001*f0);
disp(componentValueString);
componentValueString = sprintf('H = %2.2f', H);
disp(componentValueString);
componentValueString = sprintf('Q = %2.2f', Q);
disp(componentValueString);
componentValueString = sprintf('C1 = %3.3f [nF]', 1E+9*C1);
disp(componentValueString);
componentValueString = sprintf('C2 = %3.3f [nF]', 1E+9*C2);
disp(componentValueString);
componentValueString = sprintf('R1 = %3.3f [kOhms]', 0.001*R1);
disp(componentValueString);
componentValueString = sprintf('R2 = %3.3f [kOhms]', 0.001*R2);
disp(componentValueString);
componentValueString = sprintf('R3 = %3.3f [kOhms]', 0.001*R3);
disp(componentValueString);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% TRANSFER FUNCTIONS %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = linspace(1,2*f0,1E+4);
s = 1i*2*pi*f;

% Multiple Feedback BPF Transfer Function
H_MFB_BPF2 = -s./(R1*C2) ./ ...
    (s.^2+s.*(C2+C1)/(R3*C2*C1)+(R1+R2)/(R1*R2*R3*C2*C1));

% Optionally normalizing the transfer function
H_MFB_BPF2 = H_MFB_BPF2/(H*Q);

% Determing the peak value of the pass band
[~, KfPeak] = max(abs(H_MFB_BPF2));
fPeak = f(KfPeak);

% Equivalent MFB BBF Transfer Function
%H_MFB_BPF2 = -H*omega*s ./ (s.^2+2*zeta*omega.*s+omega^2);

% Determining critical points of transfer function
% sMax1 = (-2*zeta*w0-w0*sqrt(4*zeta^2-3)) / 3;
% sMax2 = (-2*zeta*w0+w0*sqrt(4*zeta^2-3)) / 3;
% 
% H_MFB_BPF2_Max1 = abs(-H*w0*sMax1 / (sMax1^2+2*zeta*w0*sMax1+w0^2));
% H_MFB_BPF2_Max2 = abs(-H*w0*sMax2 / (sMax2^2+2*zeta*w0*sMax2+w0^2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% PLOTTING TRANSFER FUNCTION %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plotting
figure(1);
    subplot(2,1,1);
        plot(0.001*f,abs(H_MFB_BPF2));% 2nd Order
        hold on;
        plot(0.001*f,abs(H_MFB_BPF2.^2));%  4th Order
        plot(0.001*f,abs(H_MFB_BPF2.^3));%  6th Order
        plot(0.001*f,abs(H_MFB_BPF2.^4));%  8th Order
        plot(0.001*f,abs(H_MFB_BPF2.^5));% 10th Order
        grid on;
        ylim([0,1]);
        xlabel('Frequency [kHz]');
        ylabel('|| H(j\omega) ||');
        legend({'2nd','4th','6th','8th','10th'});
        titleString = sprintf('NORMALIZED Bandpass Filter Frequency Response\nH = %3.1f\t\t\tQ = %3.1f\t\t\tf_{peak} = %3.3f [kHz]', H,Q,0.001*fPeak);
        title(titleString);
        hold off;
	subplot(2,1,2);
        plot(0.001*f,(180/pi)*angle(H_MFB_BPF2));% 2nd Order
        hold on;
        plot(0.001*f,(180/pi)*angle(H_MFB_BPF2.^2));%  4th Order
        plot(0.001*f,(180/pi)*angle(H_MFB_BPF2.^3));%  6th Order
        plot(0.001*f,(180/pi)*angle(H_MFB_BPF2.^4));%  7th Order
        plot(0.001*f,(180/pi)*angle(H_MFB_BPF2.^5));% 10th Order
        grid on;
        ylim([-180,180]);
        set(gca,'YTick',[-180,-135,-90,-45,0,45,90,135,180])
        set(gca,'YTickLabel',[-180,-135,-90,-45,0,45,90,135,180]);
        xlabel('Frequency [kHz]');
        ylabel('\Theta H(j\omega) [Deg]');
        legend({'2nd','4th','6th','8th','10th'});
        hold off;
        