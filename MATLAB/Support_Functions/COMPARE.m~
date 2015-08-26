function closest = COMPARE(testValue,testArray)
    % ***INCOMPLETE***
    %
    % This function compares a test value with values in an array.
    % The value in the array that is the closest to the test value
    % is returned.

    closest = testArray(1);
    diff1 = testArray(1) - testValue;

    for i=2:length(testArray);
        diff2 = testArray(i) - testValue;

        if ( abs(diff2) < abs(diff1) )
            diff1 = diff2;
            closest = testArray(i);
        end
    end
end

