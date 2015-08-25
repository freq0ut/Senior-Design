function A = SIMPSON(X,Y)
    % ***COMPLETED***
    % This function uses numerical integration using Simpsons Rule.
    % dx IS ASSUMED TO BE CONSTANT!
    
    if ( length(X) ~= length(Y) );
        error('X and Y are not the same length');
    elseif ( mod(length(X),2) ~= 0 )
        error('X and Y are not of even length');
    end
    
    % Initialization of variables
    dx = X(2)-X(1);
    A = Y(1) + Y(end);
    
    for i=2:2:length(X)-1;
        A = A + 4*Y(i);
    end
    
    for i=3:2:length(X)-1;
        A = A + 2*Y(i);
    end
    
    A = (dx/3)*A;   
end
    