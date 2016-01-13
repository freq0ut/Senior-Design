function [DEC_NUM] = BIN_2_DEC(BIN_NUM)
%   [DEC_NUM] = BIN_2_DEC(BIN_NUM)
%
%   BIN_2_DEC converts a binary number (2's Comp) represented as an array
%   of integers to a decimal number.
    
    N = length(BIN_NUM);

    DEC_NUM = 0;

    for i=0:N-2;
        DEC_NUM = DEC_NUM + BIN_NUM(N-i) * 2^i;
    end

    if ( BIN_NUM(1) == 1 )
        DEC_NUM = DEC_NUM - 2^(N-1);
    end
    
end
