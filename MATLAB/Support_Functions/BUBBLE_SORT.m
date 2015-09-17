function [Y,X] = BUBBLE_SORT(Y,X)
    % ***COMPLETED***
    % This function sorts Y from smallest to largest
    % while keeping the pairing of the X's intact.
        error('X and Y are not the same length!')
    end

    lastEntry = length(Y);

    while ( lastEntry > 0 )

        for i=1:lastEntry-1;

            if ( Y(i) > Y(i+1) )
                temp = Y(i+1);
                Y(i+1) = Y(i);
                Y(i) = temp;

                temp = X(i+1);
                X(i+1) = X(i);
                X(i) = temp;
            end
        end

	    lastEntry = lastEntry - 1;
    end
end
