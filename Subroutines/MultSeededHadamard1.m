function [Z] = MultSeededHadamard1(X, Jup, Jdown, numBlock, Mblock, Nblock, rpUp, rpMiddle, rpDown)
% Nbloc, size of blocks, must be a power of 2
% rp*, random permutation vector for the blocks

if (isrow(X) ); X = X'; end
Z = zeros(sum(Mblock), 1);

% first block line
Y2 = hadamards(X(1 : Nblock) );
Y3 = Jup * hadamards(X(1 + Nblock : 2 * Nblock) );
Z(1 : Mblock(1) ) = Y2(rpMiddle(1 : Mblock(1) ) ) + Y3(rpUp(1 : Mblock(1) ) );
lastZ = Mblock(1);

% middle block lines
aX = 1;
for block = 2 : numBlock - 1
    Y1 = Jdown * hadamards(X(aX : aX + Nblock - 1) );
    Y2 = hadamards(X(aX + Nblock : aX + 2 * Nblock - 1) );
    Y3 = Jup * hadamards(X(aX + 2 * Nblock : aX + 3 * Nblock - 1) );      
    Z(lastZ + 1 : lastZ + Mblock(block) ) = Y1(rpDown((block - 2) * Nblock + 1 : (block - 2) * Nblock + Mblock(block) ) ) + Y2(rpMiddle((block - 1) * Nblock + 1 : (block - 1) * Nblock + Mblock(block) ) ) + Y3(rpUp((block - 1) * Nblock + 1 : (block - 1) * Nblock + Mblock(block) ) );
    aX = aX + Nblock;
    lastZ = lastZ + Mblock(block);
end

% last block line
aX = (numBlock - 2) * Nblock + 1;
Y1 = Jdown * hadamards(X(aX : aX + Nblock - 1) );
Y2 = hadamards(X(aX + Nblock : numBlock * Nblock) );
Z(lastZ + 1 : lastZ + Mblock(numBlock) ) = Y1(rpDown((numBlock - 2) * Nblock + 1 : (numBlock - 2) * Nblock + Mblock(numBlock) ) ) + Y2(rpMiddle((numBlock - 1) * Nblock + 1 : (numBlock - 1) * Nblock + Mblock(numBlock) ) );

if (isrow(X) ); Z = Z'; end

end

