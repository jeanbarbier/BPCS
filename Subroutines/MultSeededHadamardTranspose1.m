function [Z] = MultSeededHadamardTranspose1(X, Jup, Jdown, numBlock, Mblock, Nblock, rpUp, rpMiddle, rpDown)
% Nbloc, size of blocks, must be a power of 2
% rp*, random permutation vector for the blocks

if (isrow(X) ); X = X'; X = [X; zeros(Nblock, 1)];
else X = [X; zeros(Nblock, 1)]; end
Z = zeros(numBlock * Nblock, 1);
zero = zeros(Nblock, 1);

% first block line
XX = zero;
S = X(1 : Nblock);
XX(rpMiddle(1 : Mblock(1) ), 1) = S(1 : Mblock(1) );
Y2 = hadamards(XX);

XX = zero;
S = X(1 + Mblock(1) : Mblock(1) + Nblock);
XX(rpDown(1 : Mblock(2) ), 1) = S(1 : Mblock(2) );
Y3 = Jdown * hadamards(XX);

Z(1 : Nblock) = Y2 + Y3;

% second block line
XX = zero;
S = X(1 : Nblock);
XX(rpUp(1 : Mblock(1) ), 1) = S(1 : Mblock(1) );
Y1 = Jup * hadamards(XX);

XX = zero;
S = X(1 + Mblock(1) : Mblock(1) + Nblock);
XX(rpMiddle(Nblock + 1 : Nblock + Mblock(2) ), 1) = S(1 : Mblock(2) );
Y2 = hadamards(XX);

XX = zero;
S = X(1 + Mblock(1) + Mblock(2) : Mblock(1) + Mblock(2) + Nblock);
XX(rpDown(Nblock + 1 : Nblock + Mblock(2) ), 1) = S(1 : Mblock(2) );
Y3 = Jdown * hadamards(XX);

Z(Nblock + 1 : 2 * Nblock) = Y1 + Y2 + Y3;

% middle block lines
aX = Mblock(1);
for block = 3 : numBlock - 1
    XX = zero;
    S = X(aX + 1 : aX + Nblock);
    XX(rpUp((block - 2) * Nblock + 1 : (block - 2) * Nblock + Mblock(block - 1)), 1) = S(1 : Mblock(block - 1) );
    Y1 = Jup * hadamards(XX);
    
    XX = zero;
    S = X(aX + Mblock(block - 1) + 1 : aX + Mblock(block - 1) + Nblock);
    XX(rpMiddle((block - 1) * Nblock + 1 : (block - 1) * Nblock + Mblock(block) ), 1) = S(1 : Mblock(block) );
    Y2 = hadamards(XX);
    
    XX = zero;
    S = X(aX + Mblock(block - 1) + Mblock(block) + 1 : aX + Mblock(block - 1) + Mblock(block) + Nblock);
    XX(rpDown((block - 1) * Nblock + 1 : (block - 1) * Nblock + Mblock(block + 1) ), 1) = S(1 : Mblock(block + 1) );
    Y3 = Jdown * hadamards(XX);
    
    Z((block - 1) * Nblock + 1 : block * Nblock) = Y1 + Y2 + Y3;
    
    aX = aX + Mblock(block - 1);
end

% last block line
XX = zero;
S = X(aX + 1 : aX + Nblock);
XX(rpUp((numBlock - 2) * Nblock + 1 : (numBlock - 2) * Nblock + Mblock(end - 1)), 1) = S(1 : Mblock(end - 1) );
Y1 = Jup * hadamards(XX);

XX = zero;
S = X(aX + Mblock(end - 1) + 1 : aX + Mblock(end - 1) + Nblock);
XX(rpMiddle((numBlock - 1) * Nblock + 1 : (numBlock - 1) * Nblock + Mblock(end)), 1) = S(1 : Mblock(end) );
Y2 = hadamards(XX);

Z((numBlock - 1) * Nblock + 1 : numBlock * Nblock) = Y1 + Y2;

if (isrow(X) ); Z = Z'; end

end

