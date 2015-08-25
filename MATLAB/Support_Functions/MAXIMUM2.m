function [yMax, xMax] = MAXIMUM2(X,Y,x1,x2)
    % ***COMPLETED***
    % This function calculates the maximum Y value (yMax) and its
    % corresponding X location (xMax).
    
    if length(X) ~= length(Y);
        error('X and Y are not the same length');
    end
    
    yMax = 0;
    
    for x=x1:x2;
       yTest = Y(x);
       
       if (yTest>yMax)
           yMax = yTest;
           xMax = x;
       end
    end
end
