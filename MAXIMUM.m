function [yMax, xMax] = MAXIMUM(X,Y)
    % This function calculates the maximum Y value (yMax) and its
    % corresponding X location (xMax).
    
    if length(X) ~= length(Y);
        error('X and Y are not the same length');
    end
    
    yMax = 0;
    
    for x=1:length(X);
       yTest = Y(x);
       
       if (yTest>yMax)
           yMax = yTest;
           xMax = x;
       end
    end
end
