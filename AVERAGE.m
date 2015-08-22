function yAvg = AVERAGE(X,Y)
    % This function calculates average y-value.
    % dx is ASSUMED TO BE CONSTANT!
    
    if length(X) ~= length(Y);
        error('X and Y are not the same length');
    end
    
    A = TRAPZ(X,Y);
    yAvg = A / ( X(end)-X(1) );
end
