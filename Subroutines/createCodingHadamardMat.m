function [G] = createCodingHadamardMat(Mblock, Nblock, numBlockL, numBlockC, rp, J, noBlockError)

G = zeros(sum(Mblock), numBlockC * Nblock);
h = hadamard(Nblock);
Y = 1;

for l = 1 : numBlockL
    
    X = 1;
    for c = 1 : numBlockC
        
        if (isempty(rp{l, c} ) == 0)
            
            hh = h(rp{l, c}, :);
            hh(1 : Mblock(l), :);
            G(Y : Y + Mblock(l) - 1, X : X + Nblock - 1) = J(l,c) * hh(1 : Mblock(l), :);
            if (max(size(noBlockError) ) > 0);
                GG = G(Y : Y + Mblock(l) - 1, X : X + Nblock - 1);
                GG(noBlockError{l}, :) = -GG(noBlockError{l}, :);
                G(Y : Y + Mblock(l) - 1, X : X + Nblock - 1) = GG;
            end
            
        end
        
        X = X + Nblock;
    end
    
    Y = Y + Mblock(l);
end

end

