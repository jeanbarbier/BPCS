function [rperm, flipedSigns] = createRandomLinesAndSignsPermutationForOperators(numBlockC, numBlockL, J, Mblock, Nblock)
% create the random permutations of the lines for the fourier or hadamard
% operators avoiding repetitions

% selection of the lines of the matrices that are sign-switched
for l = 1 : numBlockL
    flipedSigns{l} = randperm(Mblock(l), floor(Mblock(l) ./ 2) );
end

% selection of the operator modes for each column of the matrices
for c = 1 : numBlockC
    fullPerm = randperm(Nblock);
    i1 = find(fullPerm == 1);
    % No first mode
    fullPerm(i1) = fullPerm(i1 + 1);
    rpC = fullPerm;
    
    u = 1;
    for l = 1 : numBlockL
        if (J(l, c) ~= 0);
            if (u >= max(size(rpC) ) - Mblock(1) ); u = 1; rpC = fullPerm(randperm(Nblock) ); end
            rperm{l, c} = rpC(u : u + Mblock(l) - 1); u = u + Mblock(l);
        else rperm{l, c} = []; end
    end
    
end

end

