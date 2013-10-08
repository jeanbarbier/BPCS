function [Z] = MultSeededFourierInverse(X, J, numBlockL, numBlockC, Mblock, Nblock, rp, noBlockError)
% Nbloc, size of blocks, must be a power of 2
% rp*, random permutation vector for the blocks

if (isrow(X) ); X = X.'; end
Z = zeros(numBlockC * Nblock, 1);
ZZ = Z;

lastZ = 0; zero = zeros(1, Nblock);
for c = 1 : numBlockC
    
    Y = 0; YY = 0; stop = Mblock(1); start = 1;
    for l = 1 : numBlockL
        
        if (J(l, c) ~= 0)
            XXX = zero;
            XXX(rp{l, c}(1 : Mblock(l) ) ) = X(start : stop);
            Y = Y + ifft(XXX) .* Nblock .* J(l, c);
            
            if (max(size(noBlockError) ) > 0)
                XXXX = zero;
                XXXX(rp{l, c}(noBlockError{l} ) ) = XXX(rp{l, c}(noBlockError{l} ) );
                YY = YY + ifft(XXXX) .* Nblock .* J(l, c);
            end
            
        end
        
        if (l < numBlockL); start = stop + 1; stop = stop + Mblock(l + 1);  end
        
    end
    
    Z(lastZ + 1 : lastZ + Nblock) = Y;
    
    if (max(size(noBlockError) ) > 0) ZZ(lastZ + 1 : lastZ + Nblock) = YY; end
    
    lastZ = lastZ + Nblock;
end

if (max(size(noBlockError) ) > 0); Z = Z - 2 * ZZ; end

if (isrow(X) ); Z = Z.'; end

end