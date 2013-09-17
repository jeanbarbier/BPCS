function [MSE] = MSEbyBloc(X, S, numBlockC, Nblock)

MSE = zeros(1, numBlockC);
for l = 1 : numBlockC
    MSE(l) = mean((X((l - 1) * Nblock + 1: l * Nblock) - S((l - 1) * Nblock + 1: l * Nblock) ).^2);
end

end