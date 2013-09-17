function [Z] = MultSeededHadamardTransposeSquarred2(X, J, numBlockL, numBlockC, Mblock, Nblock)
% Nbloc, size of blocks, must be a power of 2
% rp*, random permutation vector for the blocks

if (isrow(X) ); X = X'; X = [X; zeros(Nblock, 1)];
else X = [X; zeros(Nblock, 1)]; end
Z = zeros(numBlockC * Nblock, 1);

lastZ = 1;
for c = 1 : numBlockC
    
    if (c > 2); u = sum(Mblock(1 : c - 2) ) + 1;
    else u = 1; end;
    
    Y = 0;
    for l = c - 1 : numBlockL
        
        if (c == 1); l = l + 1; l = min(l, numBlockL); end
        
        if (J(l, c) == 0); break; end
        
        Y = Y + J(l, c)^2 * sum(X(u : u + Mblock(l) - 1) ) * ones(Nblock, 1);
        u = u + Mblock(l);
    end
    
    Z(lastZ : lastZ + Nblock - 1) = Y;
    lastZ = lastZ + Nblock;
end

if (isrow(X) ); Z = Z'; end

end