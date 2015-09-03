function [XC, XC_Lags] = XCORR2(Y1,Y2,iBound)
    % ***COMPLETED***
    % This function calculates cross-correlation.
    
    if length(Y1) ~= length(Y2);
        error('Y1 and Y2 are not the same length');
    end
    
    N0 = length(Y1);
    X = 1:N0; % X array for TRAPZ to integrate properly
    XC = zeros(1,2*iBound+1);
    XC_Lags = zeros(1,2*iBound+1);
    Y2 = padarray(Y2,[0,iBound],0,'both'); % Y2 is padded with zeros on
                                           % both sides by iBound.
    Y2copy = zeros(1,N0); % Copy of Y2 that will be shifted continuously.

    for tau = -iBound:iBound;
        for shift = 1:N0;
           Y2copy(shift) = Y2(iBound-tau+shift);
        end
        
        A = TRAPZ(X,Y1.*Y2copy);
        XC(tau+iBound+1) = A;
        XC_Lags(tau+iBound+1) = tau;
    end
end

% Attempted to access Y2(1085); index out of bounds because numel(Y2)=1084.
%
% N0 = 20
% iBound = 5
%
% X = 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
%
% Y1 =               13 47  69 -13  58 68 15 93  1 74 *** 58 65  1  6 7 84 84  45 32  5 
% Y2 =     0 0 0 0 0 18 57 -99   4 -28 17 77 68 14 51 *** 42  8 -6 47 6  0  5 -19 58 13 0 0 0 0 0
% Y2copy =           18 57 -99   4 -28 17 77 68 14 51 *** 42  8 -6 47 6  0  5 -19 58 13
%
% Y2copy( tau = -5 ) =  17  77  68  14  51  42   8  -6  47   6   0   0   0   0   0
% Y2copy( tau = -4 ) = -28  17  77  68  14  51  42   8  -6  47   6   0   0   0   0
% Y2copy( tau = -3 ) =   4 -28  17  77  68  14  51  42   8  -6  47   6   0   0   0
% Y2copy( tau = -2 ) = -99   4 -28  17  77  68  14  51  42   8  -6  47   6   0   0
% Y2copy( tau = -1 ) =  57 -99   4 -28  17  77  68  14  51  42   8  -6  47   6   0
% Y2copy( tau =  0 ) =  18  57 -99   4 -28  17  77  68  14  51  42   8  -6  47   6
% Y2copy( tau =  1 ) =   0  18  57 -99   4 -28  17  77  68  14  51  42   8  -6  47
% Y2copy( tau =  2 ) =   0   0  18  57 -99   4 -28  17  77  68  14  51  42   8  -6
% Y2copy( tau =  3 ) =   0   0   0  18  57 -99   4 -28  17  77  68  14  51  42   8
% Y2copy( tau =  4 ) =   0   0   0   0  18  57 -99   4 -28  17  77  68  14  51  42
% Y2copy( tau =  5 ) =   0   0   0   0   0  18  57 -99   4 -28  17  77  68  14  51  
%
% length(XC) = 2*iBound+1
%
