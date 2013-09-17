function [J] = createSeededHadamardJ(numBlockL, numBlockC, JJ, w)

J = zeros(numBlockL, numBlockC);
ww = [0 : w];

for l = 1 : numBlockL
    
    for c = 1 : numBlockC
        if (c == l + 1);
            J(l, c) = JJ;
        elseif (sum(c == l - ww) > 0)
            J(l, c) = 1;
        end
    end
    
end

end

