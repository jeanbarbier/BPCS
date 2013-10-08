function [Z] = MultSeededHadamardTranspose2(X, J, numBlockL, numBlockC, Mblock, Nblock, rp, noBlockError)
% Nbloc, size of blocks, must be a power of 2
% rp*, random permutation vector for the blocks

if (isrow(X) ); X = X'; X = [X; zeros(10 * Nblock, 1)];
else X = [X; zeros(10 * Nblock, 1)]; end

Z = zeros(numBlockC * Nblock, 1);
ZZ = Z;
zero = zeros(Nblock, 1);

lastZ = 1;
for c = 1 : numBlockC
    
    if (c > 2); u = sum(Mblock(1 : c - 2) ) + 1;
    else u = 1; end;
    
    Y = 0; YY = 0;
    for l = c - 1 : numBlockL
        
        if (c == 1); l = l + 1; l = min(l, numBlockL); end
        
        if (J(l, c) == 0); break; end
        
        XX = zero;
        S = X(u : u + Nblock - 1);
        XX(rp{l, c}(1 : Mblock(l) ) ) = S(1 : Mblock(l) );
        
        Y = Y + J(l, c) * hadamards(XX);
        u = u + Mblock(l);
        
        if (max(size(noBlockError) ) > 0)
            XXX = zero;
            XXX(rp{l, c}(noBlockError{l} ) ) = XX(rp{l, c}(noBlockError{l} ) );
            YY = YY + J(l, c) * hadamards(XXX);
        end
        
    end
    
    Z(lastZ : lastZ + Nblock - 1) = Y;
    
    if (max(size(noBlockError) ) > 0) ZZ(lastZ : lastZ + Nblock - 1) = YY; end
    
    lastZ = lastZ + Nblock;
end

if (max(size(noBlockError) ) > 0); Z = Z - 2 * ZZ; end

if (isrow(X) ); Z = Z'; end

end