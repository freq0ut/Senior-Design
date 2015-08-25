function [pks, pkLocs] = FIND_POS_PEAKS(Y,THD,MPD,NBRS)
    % INCOMPLETE
    % This function finds the local maxima of the array Y
    % THD = Threshold
    % MPD = Minimum Peak Distance
    % NBRS = Neighbors
    
    % Finding yMax and its corresponding index i_yMax
    yMax = 0;
    for i=1:length(Y);
       yTest = Y(i);
       
       if yTest > yMax
           yMax = yTest;
           i_yMax = i;
       end
    end
    
    % Checking to see if THD is set too high
    if ( yMax < THD )
        error('No peaks were found');
    else
        pks(1) = yMax;
        pkLocs(1) = i_yMax;
        pkCount = 2;
    end

    % Looking to the right of i_yMax
    i = i_yMax+MPD;
    while i < length(Y)
        yTest = Y(i);
        
        % Left neighbors
        iL = i-NBRS:i-1;
        
        % Right neighbors
        iR = i+1:i+NBRS;
        
        if (iL(1) > 0 && iR(end) < length(Y)-1)
            leftNBRS  = Y(iL);
            rightNBRS = Y(iR);
            
            if ( yTest>THD && yTest>max(leftNBRS) && yTest>max(rightNBRS) )
                pks(pkCount) = yTest;
                pkLocs(pkCount) = i;
                i = i + MPD;
                pkCount = pkCount + 1;
            end
        end
        
        i = i+1;
    end
    
    % Looking to the left of i_yMax
end
