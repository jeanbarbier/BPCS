function [Z] = MultSeededHadamardSquarred1(X, Jup, Jdown, numBlock, Mblock, Nblock)
% Nbloc, size of blocks, must be a power of 2
% rp*, random permutation vector for the blocks

if (isrow(X) ); X = X'; end
Z = zeros(sum(Mblock), 1);

% first block line
Y2 = sum(X(1 : Nblock) ) * ones(Mblock(1), 1);
Y3 = Jup^2 * sum(X(1 + Nblock : 2 * Nblock) ) * ones(Mblock(1), 1);
Z(1 : Mblock(1) ) = Y2 + Y3;
lastZ = Mblock(1);

% middle block lines
aX = 1;
for block = 2 : numBlock - 1
    Y1 = Jdown^2 * sum(X(aX : aX + Nblock - 1) ) * ones(Mblock(block), 1);
    Y2 = sum(X(aX + Nblock : aX + 2 * Nblock - 1) ) * ones(Mblock(block), 1);
    Y3 = Jup^2 * sum(X(aX + 2 * Nblock : aX + 3 * Nblock - 1) ) * ones(Mblock(block), 1);
    Z(lastZ + 1 : lastZ + Mblock(block), 1) = Y1 + Y2 + Y3;
    aX = aX + Nblock;
    lastZ = lastZ + Mblock(block);
end

% last block line
aX = (numBlock - 2) * Nblock + 1;
Y1 = Jdown^2 * sum(X(aX : aX + Nblock - 1) ) * ones(Mblock(numBlock), 1);
Y2 = sum(X(aX + Nblock : numBlock * Nblock) ) * ones(Mblock(numBlock), 1);
Z(lastZ + 1 : lastZ + Mblock(numBlock) ) = Y1 + Y2;

if (isrow(X) ); Z = Z'; end

end

