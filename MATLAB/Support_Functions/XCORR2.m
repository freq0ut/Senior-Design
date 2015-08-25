function [XC, XC_Lags] = XCORR2(Y1,Y2,i)
    % ***INCOMPLETE***
    % This function calculates cross-correlation.
    
    if length(Y1) ~= length(Y2);
        error('Y1 and Y2 are not the same length');
    end
    
    N0 = length(Y1);
    X = 1:N0; % X array for TRAPZ to integrate properly
    XC = zeros(1,2*N0-1);
    XC_Lags = zeros(1,2*N0-1);    
    Y2 = padarray(Y2,[0,N0],0,'both'); % Y2 is padded with zeros on
                                       % both sides TRIPLING its size.
    Y2s = zeros(1,N0); % Smaller Y2 array that will be shifted continuously.

    for tau = -i:i;
        for shift = 1:N0; % Shifting the smaller Y2 array
           Y2s(shift) = Y2(shift-tau+N0);
        end
        
        A = TRAPZ(X,Y1.*Y2s);
        XC(N0+tau) = A;
        XC_Lags(N0+tau) = tau;
    end
end
