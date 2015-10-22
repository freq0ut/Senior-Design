% Author: Joshua Simmons
%
% Date: October 22, 2015
%
% Desciption: Selects the R and C values as well as plot the freq response for
%             the Multiple Feedback Bandpass Filter. 

clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% USER INPUT %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fMid = 20E+3;% Mid frequency

% THE CIRCUIT GAIN IS THE PRODUCT OF H AND Q
H = 10;% Gain at mid freq
Q = 10;% Quality factor

% Select the value of C3 to give you the desired freq response.
C3 = 7.95E-6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% CALCULATING COMPONENT VALUES %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

omega = 2*pi*fMid;

zeta = 1/(2*Q);    % Damping ratio
alpha = omega*zeta;% Normalized damping ratio

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

% Multiple Feedback BPF Transfer Function
H_MFB_BPF2 = -s.*(1/(R1*C4)) ./ ...
    (s.^2+s.*(C3+C4)/(C3*C4*R5)+(1/(R5*C3*C4))*(1/R1+1/R2));

% Equivalent MFB BBF Transfer Function
%H_MFB_BPF2 = -H*omega*s ./ (s.^2+2*zeta*omega.*s+omega^2);

% Normalized Prototypical 2nd Order BPF Transfer Function
H_PRO_BPF2 = (1/Q) * (omega^2) ./ (s.^2+2*zeta*omega.*s+omega^2);

% Optionally normalizing the transfer function
H_MFB_BPF2 = H_MFB_BPF2/(H*Q);

% Determing the peak value of the pass band
[~, KfPeak] = max(abs(H_MFB_BPF2));
fPeak = f(KfPeak);

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
        plot(0.001*f,abs(H_PRO_BPF2));   %  2nd Order Proto
        grid on;
        ylim([0,1]);
        xlabel('Frequency [kHz]');
        ylabel('|| H(j\omega) ||');
        legend({'2nd','4th','6th','8th','10th','2nd Proto'});
        titleString = sprintf('Bandpass Filter Frequency Response\nH = %3.1f\t\t\tQ = %3.1f\t\t\tf_{peak} = %3.3f [kHz]', H,Q,0.001*fPeak);
        title(titleString);
        hold off;
	subplot(2,1,2);
        plot(0.001*f,(180/pi)*angle(H_MFB_BPF2));% 2nd Order
        hold on;
        plot(0.001*f,(180/pi)*angle(H_MFB_BPF2.^2));%  4th Order
        plot(0.001*f,(180/pi)*angle(H_MFB_BPF2.^3));%  6th Order
        plot(0.001*f,(180/pi)*angle(H_MFB_BPF2.^4));%  7th Order
        plot(0.001*f,(180/pi)*angle(H_MFB_BPF2.^5));% 10th Order
        plot(0.001*f,(180/pi)*angle(H_PRO_BPF2));   %  2nd Order Proto
        grid on;
        ylim([-180,180]);
        set(gca,'YTick',[-180,-135,-90,-45,0,45,90,135,180])
        set(gca,'YTickLabel',[-180,-135,-90,-45,0,45,90,135,180]);
        xlabel('Frequency [kHz]');
        ylabel('\Theta H(j\omega) [Deg]');
        legend({'2nd','4th','6th','8th','10th','2nd Proto'});
        hold off;
        