function TOA = FIND_TOA(d,td2,td3,td4,vP)
    % This function uses the Bisector Method to find the Time-Of-Arrival.
    
    function y = f(d,td2,td3,td4,TOA,vP)
        % This is the TOA function.
        
        % Parameter initialization
        R1 = vP*(TOA);
        R2 = vP*(TOA+td2);
        R3 = vP*(TOA+td3);
        R4 = vP*(TOA+td4);
        
        y =   ( -R1.^2 + R2.^2 + R3.^2 - R4.^2 - 2*d^2 ).^2 ...
            + (  R1.^2 - R2.^2 + R3.^2 - R4.^2 + 2*d^2 ).^2 ...
            + (  R1.^2 + R2.^2 - R3.^2 - R4.^2 + 2*d^2 ).^2 ...
            - (4*d*R1).^2;
    end
    
    % Bisector Method
    TOA1 = 0; % Left  most starting value
    TOA2 = 1; % Right most starting value
    for count = 1:30;
        TOA_Test = (TOA2+TOA1)/2;
        
        if ( f(d,td2,td3,td4,TOA_Test,vP) * f(d,td2,td3,td4,TOA2,vP) < 0 )
            TOA1 = TOA_Test;
        else
            TOA2 = TOA_Test;
        end
    end
    
    TOA = TOA1;

end
