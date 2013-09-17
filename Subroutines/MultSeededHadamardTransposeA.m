function [Z] = MultSeededHadamardTransposeA(X, J, numBlockL, numBlockC, Mblock, Nblock, rp, noBlockError)
% Nbloc, size of blocks, must be a power of 2
% rp*, random permutation vector for the blocks

if (isrow(X) ); X = X'; X = [X; zeros(2 * Nblock, 1)];
else X = [X; zeros(2 * Nblock, 1)]; end

Z = zeros(numBlockC * Nblock, 1);
ZZ = Z;
zero = zeros(Nblock, 1);

lastZ = 1;
for c = 1 : numBlockC
    
    u = 1;
    Y = 0;
    YY = 0;
    for l = 1 : numBlockL
        
        if (J(l, c) ~= 0)
            
            XX = zero;
            S = X(u : u + Nblock - 1);
            XX(rp{l, c}(1 : Mblock(l) ) ) = S(1 : Mblock(l) );
            
            Y = Y + J(l, c) * hadamards(XX);
            
            if (max(size(noBlockError) ) > 0)
                XXX = zero;
                XXX(rp{l,c}(noBlockError{l} ) ) = XX(rp{l,c}(noBlockError{l} ) );
                YY = YY + J(l, c) * hadamards(XXX);
            end
            
        end
        
        u = u + Mblock(l);
        
    end
    
    Z(lastZ : lastZ + Nblock - 1) = Y;
    
    if (max(size(noBlockError) ) > 0)
        ZZ(lastZ : lastZ + Nblock - 1) = YY;
    end
    
    lastZ = lastZ + Nblock;
end

if (max(size(noBlockError) ) > 0); Z = Z - 2 * ZZ; end

if (isrow(X) ); Z = Z'; end

end