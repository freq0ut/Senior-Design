% Author: Joshua Simmons
%
% Date: October 22, 2015
%
% Desciption: Selects the R and C values for a Multiple Feedback BPF.

clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% USER INPUT %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fMid = 20E+3;% Mid frequency

% NOTE: The circuit gain at fMid is H*Q
H = 10;% Partial gain at the mid frequency
Q = 10;% Quality factor

C3 = 8E-6;% You must pick C3. You will need to change it to get the
           % frequency response you like.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% CALCULATING COMPONENT VALUES %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

omega = 2*pi*fMid;
zeta = 1/(2*Q);
alpha = omega*zeta;

k = omega*C3;
C4 = C3;
R1 = 1/(H*k);
R2 = 1/((2*Q-H))*k;
R5 = 2*Q/k;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% TRANSFER FUNCTIONS %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = linspace(1,2*fMid,1E+4);
s = 1i*2*pi*f;

% Multiple Feedback BPF
H_MFB_BPF = -s.*(1/(R1*C4)) ./ ...
    (s.^2+s.*(C3+C4)/(C3*C4*R5)+(1/(R5*C3*C4))*(1/R1+1/R2));

% Equivalently
H_PRO_BPF = -H*omega*s ./ (s.^2+alpha*omega.*s+omega^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% PLOTTING TRANSFER FUNCTION %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plotting
figure(1);
    subplot(2,1,1);
        plot(0.001*f,abs(H_MFB_BPF));
        hold on;
        plot(0.001*f,abs(H_PRO_BPF));
        grid on;
        ylim([0,H*Q]);
        xlabel('Frequency [kHz]');
        ylabel('|| H(j\omega) ||');
        titleString = sprintf('Bandpass Filter Frequency Response\nH = %3.1f\nQ = %3.1f', H,Q);
        title(titleString);
        hold off;
	subplot(2,1,2);
        plot(0.001*f,(180/pi)*angle(H_MFB_BPF));
        hold on;
        plot(0.001*f,(180/pi)*angle(H_PRO_BPF));
        grid on;
        ylim([-180,180]);
        set(gca,'YTick',[-180,-135,-90,-45,0,45,90,135,180])
        set(gca,'YTickLabel',[-180,-135,-90,-45,0,45,90,135,180]);
        xlabel('Frequency [kHz]');
        ylabel('\Theta H(j\omega) [Deg]');
        hold off;
        