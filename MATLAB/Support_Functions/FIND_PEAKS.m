function [pks, pkLocs] = FIND_PEAKS(Y,MPD,xNBRS)
    % ***INCOMPLETE***
    % This function finds local maxima in specified zones.
    %
    % MPD = Minimum Peak Distance
    %
    % xNBRS = neighbors.
    % This is a zone to the left & right of local yMaxs
    % the algorithm searches in to find other local yMaxs.
    
    X=1:length(Y);
    
    % Finding yMax and its corresponding index xMax
    [yMax, xMax] = MAXIMUM(X,Y);
    pkLocs(1,1) = xMax;
    pks(1,1) = yMax;
    pkCount = 2;

    % Looking to the right of xMax
    x1  = xMax+MPD-1;
    x2 = x1+xNBRS;
    while x2 < length(Y)-1;
        [yMax, xMaxLocal] = MAXIMUM2(X,Y,x1,x2);
        pkLocs(1,pkCount) = xMaxLocal;
        pks(1,pkCount) = yMax;
        pkCount = pkCount+1;
        x1  = x2+MPD-1;
        x2 = x1+xNBRS;
    end
       
    % Looking to the left of xMax
    x2 = xMax-MPD+1;
    x1  = x2-xNBRS;
    while x1 > 0;
        [yMax, xMaxLocal] = MAXIMUM2(X,Y,x1,x2);
        pkLocs(1,pkCount) = xMaxLocal;
        pks(1,pkCount) = yMax;
        pkCount = pkCount+1;
        x2 = x1-MPD+1;
        x1 = x2-xNBRS;
    end
    
    % Sorting results from smallest to largest index
    [pkLocs,pks] = BUBBLE_SORT(pkLocs,pks);
end
