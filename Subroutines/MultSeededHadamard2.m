function [Z] = MultSeededHadamard2(X, J, numBlockL, numBlockC, Mblock, Nblock, rp, noBlockError)
% Nbloc, size of blocks, must be a power of 2
% rp*, random permutation vector for the blocks

if (isrow(X) ); X = X'; end
Z = zeros(sum(Mblock), 1);

lastZ = 0;
for l = 1 : numBlockL
    
    Y = 0;
    for c = min(l + 1, numBlockC) : -1 : 1
        if (J(l, c) == 0); break; end
        YY = J(l, c) * hadamards(X((c - 1) * Nblock + 1 : c * Nblock) );
        Y = Y + YY(rp{l, c}(1 : Mblock(l) ) );
    end
    
    if (max(size(noBlockError) ) > 0); Y(noBlockError{l} ) = -Y(noBlockError{l} ); end
    
    Z(lastZ + 1 : lastZ + Mblock(l) ) = Y;
    lastZ = lastZ + Mblock(l);
end

if (isrow(X) ); Z = Z'; end

end