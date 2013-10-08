function [MSE] = MSEbyBloc(X, S, numBlockC, Nblock, type)

if (strcmp(type, 'AMPhadamardSeeded') )
    MSE = zeros(1, numBlockC);
    for l = 1 : numBlockC
        MSE(l) = mean((X((l - 1) * Nblock + 1: l * Nblock) - S((l - 1) * Nblock + 1: l * Nblock) ).^2);
    end
else
    MSE = zeros(2, numBlockC);
    for l = 1 : numBlockC
        MSE(1, l) = mean((real(X((l - 1) * Nblock + 1: l * Nblock) ) - real(S((l - 1) * Nblock + 1: l * Nblock) ) ).^2);
        MSE(2, l) = mean((imag(X((l - 1) * Nblock + 1: l * Nblock) ) - imag(S((l - 1) * Nblock + 1: l * Nblock) ) ).^2);
    end
end

end