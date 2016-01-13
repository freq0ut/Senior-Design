function [BIN_NUM] = HEX_2_BIN(HEX_NUM)
%   [BIN_NUM] = HEX_2_BIN (HEX_NUM)
%   
%   Converts a character array of hexadecimal values to
%   a character array of binary values.

    N_HEX = length(HEX_NUM);
    N_BIN = 4*N_HEX;
    
    BIN_NUM = zeros(1,N_BIN);
    
    for i=1:N_HEX;
        i2 = 4*(i-1)+1;
        
        if ( HEX_NUM(i) == '0' )
            BIN_NUM(i2)   = 0;
            BIN_NUM(i2+1) = 0;
            BIN_NUM(i2+2) = 0;
            BIN_NUM(i2+3) = 0;
        elseif ( HEX_NUM(i) == '1' )
            BIN_NUM(i2)   = 0;
            BIN_NUM(i2+1) = 0;
            BIN_NUM(i2+2) = 0;
            BIN_NUM(i2+3) = 1;
        elseif ( HEX_NUM(i) == '2' )
            BIN_NUM(i2)   = 0;
            BIN_NUM(i2+1) = 0;
            BIN_NUM(i2+2) = 1;
            BIN_NUM(i2+3) = 0;
        elseif ( HEX_NUM(i) == '3' )
            BIN_NUM(i2)   = 0;
            BIN_NUM(i2+1) = 0;
            BIN_NUM(i2+2) = 1;
            BIN_NUM(i2+3) = 1;
        elseif ( HEX_NUM(i) == '4' )
            BIN_NUM(i2)   = 0;
            BIN_NUM(i2+1) = 1;
            BIN_NUM(i2+2) = 0;
            BIN_NUM(i2+3) = 0;
        elseif ( HEX_NUM(i) == '5' )
            BIN_NUM(i2)   = 0;
            BIN_NUM(i2+1) = 1;
            BIN_NUM(i2+2) = 0;
            BIN_NUM(i2+3) = 1;
        elseif ( HEX_NUM(i) == '6' )
            BIN_NUM(i2)   = 0;
            BIN_NUM(i2+1) = 1;
            BIN_NUM(i2+2) = 1;
            BIN_NUM(i2+3) = 0;
        elseif ( HEX_NUM(i) == '7' )
            BIN_NUM(i2)   = 0;
            BIN_NUM(i2+1) = 1;
            BIN_NUM(i2+2) = 1;
            BIN_NUM(i2+3) = 1;
        elseif ( HEX_NUM(i) == '8' )
            BIN_NUM(i2)   = 1;
            BIN_NUM(i2+1) = 0;
            BIN_NUM(i2+2) = 0;
            BIN_NUM(i2+3) = 0;
        elseif ( HEX_NUM(i) == '9' )
            BIN_NUM(i2)   = 1;
            BIN_NUM(i2+1) = 0;
            BIN_NUM(i2+2) = 0;
            BIN_NUM(i2+3) = 1;
        elseif ( HEX_NUM(i) == 'A' )
            BIN_NUM(i2)   = 1;
            BIN_NUM(i2+1) = 0;
            BIN_NUM(i2+2) = 1;
            BIN_NUM(i2+3) = 0;
        elseif ( HEX_NUM(i) == 'B' )
            BIN_NUM(i2)   = 1;
            BIN_NUM(i2+1) = 0;
            BIN_NUM(i2+2) = 1;
            BIN_NUM(i2+3) = 1;
        elseif ( HEX_NUM(i) == 'C' )
            BIN_NUM(i2)   = 1;
            BIN_NUM(i2+1) = 1;
            BIN_NUM(i2+2) = 0;
            BIN_NUM(i2+3) = 0;
        elseif ( HEX_NUM(i) == 'D' )
            BIN_NUM(i2)   = 1;
            BIN_NUM(i2+1) = 1;
            BIN_NUM(i2+2) = 0;
            BIN_NUM(i2+3) = 1;
        elseif ( HEX_NUM(i) == 'E' )
            BIN_NUM(i2)   = 1;
            BIN_NUM(i2+1) = 1;
            BIN_NUM(i2+2) = 1;
            BIN_NUM(i2+3) = 0;
        elseif ( HEX_NUM(i) == 'F' )
            BIN_NUM(i2)   = 1;
            BIN_NUM(i2+1) = 1;
            BIN_NUM(i2+2) = 1;
            BIN_NUM(i2+3) = 1;
        end
    end
end
