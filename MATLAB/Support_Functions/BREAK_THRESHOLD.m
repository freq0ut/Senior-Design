function iBreak = BREAK_THRESHOLD(Y,THD,DIR)
    % ***COMPLETED***
    %
    % FOR AC SIGNALS WITH ZERO DC OFFSET!
    %
    % This function evaluates an array Y from left-to-right or
    % right-to-left. The 1st Y that breaks a positive or 
    % negative threshold has its index returned.
    %
    % THD = Threshold
    % DIR = Direction (LR||RL)
    
    iBreak = -1;
    
    if DIR == 'LR'

        i = 1;
        while ( i < length(Y) + 1)
            yTest = abs(Y(i));

            if ( yTest > THD )
                iBreak = i;
                break;
            end

            i = i+1;
        end

    elseif DIR == 'RL'

        i = length(Y);
        while i > 0
            yTest = abs(Y(i));

            if ( yTest > THD )
                iBreak = i;
                break;
            end

            i = i-1;
        end
    else

        error('Input for direction should be LR or RL.');
    end

%     if ( iBreak == -1 )
%         disp('Threshold wasnt broken.')
%     end

end

