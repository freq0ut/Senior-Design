function [XC, XC_Lags] = XCORR3(Y1,Y2,tauStart,tauStop)
    % ***COMPLETED***
    % This function calculates cross-correlation.
    
    N = length(Y1);
    XC = zeros(1,2*N-1);
    XC_Lags = zeros(1,2*N-1);
    
    for tau = tauStart:tauStop;
       
       XC_Lags(N+tau) = tau;
       
       i = 1;
       
       if (tau<=0)
           while (-tau+i<=N)
               XC(N+tau) = XC(N+tau) + Y1(i) * Y2(-tau+i);
               i = i + 1;
           end
           
       else
           while (tau+i<=N)
               XC(N+tau) = XC(N+tau) + Y1(tau+i) * Y2(i);
               i = i + 1;
           end
           
       end
       
    end
    