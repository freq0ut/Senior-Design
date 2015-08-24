function [yMax, xMax] = FAST_MAXIMUM(X,Y,i1,i2)
    % This function calculates the maximum Y value (yMax) and its
    % corresponding X location (xMax).
    
    if length(X) ~= length(Y);
        error('X and Y are not the same length');
    end
    
    yMax = 0;
    
    for i=i1:i2;
       yTest = Y(i);
       
       if (yTest>yMax)
           yMax = yTest;
           xMax = i;
       end
    end
end
