close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% VARIABLE INITIALIZATIONS %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

chan2TailLength = 100;

SNR = 20;

fPinger = 40.0E+3;% Pinger frequency

N = 1024;% Samples to collect
fS = 900.0E+3;% Sampling Frequency
tS = 1/fS;% Sampling Period

t=0:tS:(N-1)*tS;% Time Array

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% CONSTRUCTING PINGER BURST %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tBurst=0:tS:1.3E-3;
burst = sin(2*pi*fPinger*tBurst);
burst = awgn(burst,SNR);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% CONSTRUCTING CHAN1 %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

chan1 = zeros(1,N);

chan1(1,1:N/2) = burst(1,N/2:N-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% CONSTRUCTING CHAN2 %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

chan2 = zeros(1,N);

chan2(1,1:chan2TailLength) = burst(1,N-chan2TailLength:N-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% FINDING TIME DELAY %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t1_Crude = tS*BREAK_THRESHOLD(chan1,0.5,'RL');
t2_Crude = tS*BREAK_THRESHOLD(chan2,0.5,'RL');

tD2_Crude = t1_Crude-t2_Crude;

[XC, XC_Lags] = xcorr(chan1,chan2,'coeff');
XC_Max = max(XC);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tBound = 0.05;% In milliseconds

stringA = sprintf('t_{1,Crude} = %f [ms]', 1E3*t1_Crude);
stringB = sprintf('t_{2,Crude} = %f [ms]', 1E3*t2_Crude);
stringC = sprintf('tD_{21.Crude} = %f [ms]', 1E3*tD2_Crude);

figure(1);
    subplot(2,2,1);
        stem(1E3*t,chan1,'-b');
        xlabel('Time [ms]');
        legend(stringA);
    subplot(2,2,3);
        stem(1E3*t,chan2,'-r');
        xlabel('Time [ms]');
        legend(stringB);
    subplot(2,2,2);
        stem(1E3*tS*XC_Lags,XC,'-m');
        hold on;
        plot(1E3*tD2_Crude,XC_Max,'g.','MarkerSize',20);
        xlabel('Time [ms]');
        xlim([1E3*tD2_Crude-tBound, 1E3*tD2_Crude+tBound]);
        legend(stringC);
        hold off;
    subplot(2,2,4);
        stem(1E3*tS*XC_Lags,XC,'-m');
        hold on;
        plot(1E3*tD2_Crude,XC_Max,'g.','MarkerSize',20);
        xlabel('Time [ms]');
        legend(stringC);
        hold off;
        