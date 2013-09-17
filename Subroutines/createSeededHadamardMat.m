function [G, J] = createSeededHadamardMat(Mblock, Nblock, numBlockL, numBlockC, JJ, rp, w, noBlockError)

G = zeros(sum(Mblock), numBlockC * Nblock);
h = hadamard(Nblock);
J = zeros(numBlockL, numBlockC);
Y = 1;
ww = [0 : w];

for l = 1 : numBlockL
    
    X = 1;
    for c = 1 : numBlockC
        if (c == l + 1);
            hh = JJ * h(rp{l, c}, :);
            G(Y : Y + Mblock(l) - 1, X : X + Nblock - 1) = hh(1 : Mblock(l), :);
            if (max(size(noBlockError) ) > 0);
                GG = G(Y : Y + Mblock(l) - 1, X : X + Nblock - 1);
                GG(noBlockError{l}, :) = -GG(noBlockError{l}, :);
                G(Y : Y + Mblock(l) - 1, X : X + Nblock - 1) = GG;
            end
            
            J(l, c) = JJ;
        elseif (sum(c == l - ww) > 0)
            l = min(l, numBlockL);
            c = min(l, numBlockC);
            hh = h(rp{l, c}, :);
            G(Y : Y + Mblock(l) - 1, X : X + Nblock - 1) = hh(1 : Mblock(l), :);
            if (max(size(noBlockError) ) > 0);
                GG = G(Y : Y + Mblock(l) - 1, X : X + Nblock - 1);
                GG(noBlockError{l}, :) = -GG(noBlockError{l}, :);
                G(Y : Y + Mblock(l) - 1, X : X + Nblock - 1) = GG;
            end
            J(l, c) = 1;
        end
        
        X = X + Nblock;
    end
    
    Y = Y + Mblock(l);
end

end

