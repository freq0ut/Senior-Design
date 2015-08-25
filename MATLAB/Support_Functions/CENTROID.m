function [XBar, YBar, A] = CENTROID(X,Y)
    % ***COMPLETED***
    % This function calculates the centroid of a lamina using the
    % Trapezoidal Rule.
    
    if length(X) ~= length(Y);
        error('X and Y are not the same length');
    end
    
    % Initialization of variables
    A  = 0;
    Sx = 0;
    Sy = 0;
    dx = X(2)-X(1); % Assumes dx is consant!
        
    for k=1:length(X)-1;
        if ( Y(k)*Y(k+1) > 0 )
            % Areas [units^2]
            a = (dx/2) * ( Y(k)+Y(k+1) );
            A = A + a;

            if ( a ~=0)
                % 1st Moment-of-Area about x-axis [units^3]
                sx = (dx*dx/6) * ( Y(k)+2*Y(k+1) );
                xbar = sx/a + X(k);
                Sx = Sx + xbar*a;

                % 1st Moment-of-Area about y-axis [units^3]
                sy = (dx/6) * ( Y(k)^2+Y(k)*Y(k+1)+Y(k+1)^2 );
                ybar = sy/a + 0;
                Sy = Sy + ybar*a;
            end
        else
            % Triangle bottom legs
            m = (Y(k+1)-Y(k))/dx;
            
            if (m ~= 0)
                xL = (m*0-Y(k))/m;
                xR = dx - xL;
            
                % Areas [units^2]
                a1 = (xL/2) * Y(k);
                a2 = (xR/2) * Y(k+1);
                A = A + a1 + a2;

                if ( a1 ~= 0 && a2 ~= 0)
                    % 1st Moment-of-Area about x-axis [units^3]
                    sx1 = (xL^2/6)*Y(k);
                    sx2 = (xR^2/6)*Y(k+1);
                    xbar = sx1/a1 + X(k) + sx2/a2 + X(k) + xL;
                    Sx = Sx + xbar*(a1+a2);

                    % 1st Moment-of-Area about y-axis [units^3]
                    sy1 = (xL/6)*Y(k)^2;
                    sy2 = (xR/6)*Y(k+1)^2;
                    ybar = sy1/a1 + 0 + sy2/a2 + 0;
                    Sy = Sy + ybar*(a1+a2);
                end
            end
        end
    end
        
    if (abs(A) > 0)
        XBar = Sx/A;
        YBar = Sy/A;
    else
        XBar = 0;
        YBar = 0;
    end
end