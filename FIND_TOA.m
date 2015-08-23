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

%     % Change of Sign Method
% roots = zeros(1,4);
%     i = 1;
%     for TOA1 = 0:1e-3:TOA_Upper;
%         TOA2 = TOA1+1e-3;
%         if ( f(d,td2,td3,td4,TOA1,vP) * f(d,td2,td3,td4,TOA2,vP) < 0 )
%             roots(i)   = TOA1;
%             roots(i+1) = TOA2;
%             i = i+2;
%         end
%     end
%     
%     % Getting ready for Bisector Method
%     if ( roots(4) > roots (2) )
%         TOA1 = roots(3);
%         TOA2 = roots(4);
%     elseif ( roots(4) < roots (2) )
%         TOA1 = roots(1);
%         TOA2 = roots(2);
%     else
%         TOA = -1;
%     end
    
%     TOA_Test = 0:1e-3:TOA_Upper;
%     Y = f(d,td2,td3,td4,TOA_Test,vP);
%         figure(99)
%         plot(TOA_Test,Y,'-r');
%         hold on;
%         line([TOA_Test(1),TOA_Test(end)],[0,0],'Color',[0,0,1]);
%         plot(TOA,f(d,td2,td3,td4,TOA,vP),'g.','MarkerSize',20);
%         string123 = sprintf('%3.3e', TOA);
%         legend(string123);
%         hold off;
end
