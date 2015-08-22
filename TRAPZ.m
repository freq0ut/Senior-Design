function A = TRAPZ(X,Y)
    % This function uses numerical integration using the trapezoidal rule.
    % dx IS ASSUMED TO BE CONSTANT!
    
    if length(X) ~= length(Y);
        error('X and Y are not the same length');
    end
    
    % Initialization of variables
    A  = 0;
    dx = X(2)-X(1); % Assumes dx is constant!
        
    for i=1:length(X)-1;
        
        % Trapezoids
        if ( Y(i)*Y(i+1) > 0 )
            a = (dx/2) * ( Y(i)+Y(i+1) );
            A = A + a;
            
        % Two triangles
        else
            m = (Y(i+1)-Y(i))/dx; % Slope of line
            
            if (m ~= 0) % NaN protection
                xL = (m*0-Y(i))/m; % Left  triangle bottom leg
                xR = dx - xL;      % Right triangle bottom leg
            
                % Areas [units^2]
                a  = (xL/2) * Y(i);
                a2 = (xR/2) * Y(i+1);
                A = A + a + a2;
            end
        end
    end
end
