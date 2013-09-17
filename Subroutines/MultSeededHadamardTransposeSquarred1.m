function [Z] = MultSeededHadamardTransposeSquarred1(X, Jup, Jdown, numBlock, Mblock, Nblock)
% Nbloc, size of blocks, must be a power of 2
% rp*, random permutation vector for the blocks

if (isrow(X) ); X = X'; X = [X; zeros(Nblock, 1)];
else X = [X; zeros(Nblock, 1)]; end
Z = zeros(numBlock * Nblock, 1);

% first block line
Y2 = sum(X(1 : Nblock) );
Y3 = Jdown^2 * sum(X(1 + Mblock(1) : Mblock(1) + Nblock) );
Z(1 : Nblock) = Y2 + Y3;

% second block line
Y1 = Jup^2 * sum(X(1 : Nblock) );
Y2 = sum(X(1 + Mblock(1) : Mblock(1) + Nblock) );
Y3 = Jdown^2 * sum(X(1 + Mblock(1) + Mblock(2) : Mblock(1) + Mblock(2) + Nblock) );
Z(Nblock + 1 : 2 * Nblock) = Y1 + Y2 + Y3;

% middle block lines
aX = Mblock(1);
for block = 3 : numBlock - 1
    Y1 = Jup^2 * sum(X(aX + 1 : aX + Nblock) );
    Y2 = sum(X(aX + Mblock(block - 1) + 1 : aX + Mblock(block - 1) + Nblock) );
    Y3 = Jdown^2 * sum(X(aX + Mblock(block - 1) + Mblock(block) + 1 : aX + Mblock(block - 1) + Mblock(block) + Nblock) );
    Z((block - 1) * Nblock + 1 : block * Nblock) = Y1 + Y2 + Y3;
    aX = aX + Mblock(block - 1);
end

% last block line
Y1 = Jup^2 * sum(X(aX + 1 : aX + Nblock));
Y2 = sum(X(aX + Mblock(end - 1) + 1 : aX + Mblock(end - 1) + Nblock) );
Z((numBlock - 1) * Nblock + 1 : numBlock * Nblock) = Y1 + Y2;

if (isrow(X) ); Z = Z'; end

end

