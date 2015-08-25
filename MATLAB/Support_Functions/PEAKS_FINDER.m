function [pks, pkLocs] = PEAKS_FINDER(Y,xNBRS)
    % ***INCOMPLETE***
    % This function finds local maxima in specified zones.
    % xNBRS = neighbors. This is a zone to the right and to the left
    % of local yMaxs that local maxes are computed.
    
    X=1:length(Y);
    
    % Finding yMax and its corresponding index xMax
    [yMax, xMax] = MAXIMUM(X,Y);
    pks(1,1) = xMax;
    pks(1,2) = yMax;
    pkCount = 3;

    % Looking to the right of xMax
    x  = xMax+1;
    x2 = x+xNBRS;
    while x2 < length(Y)-1;
        [yMax, xMaxLocal] = FAST_MAXIMUM(X,Y,x,x2);
        pks(1,pkCount)   = xMaxLocal;
        pks(1,pkCount+1) = yMax;
        pkCount = pkCount+2;
        x  = x2;
        x2 = x+xNBRS;
        disp('Test');
    end
       
    % Looking to the left of xMax
    x2 = xMax-1;
    x  = x2-xNBRS;
    while x > 0;
        [yMax, xMaxLocal] = FAST_MAXIMUM(X,Y,x,x2);
        pks(1,pkCount)   = xMaxLocal;
        pks(1,pkCount+1) = yMax;
        pkCount = pkCount+2;
        x2 = x;
        x  = x2-xNBRS;
    end
    
    %FAST_BUBBLE_SORT(pks,pkLocs)
end
