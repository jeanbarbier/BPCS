clear

N = 200;
M = 100;
numBlockC = 5;
numBlockL = 4;
Nblock = N ./ numBlockC;

SN = (randn(1, N) + 1i .* randn(1, N) )';
SM = (randn(1, M) + 1i .* randn(1, M) )';


for l = 1 : numBlockL
    Mblock(l) = M ./ numBlockL;
    noBlockError{l} = randperm(Mblock(l), floor(Mblock(l) ./ 2) );
%     noBlockError{l} = [];
    for c = 1 : numBlockC
        rp{l, c} = randperm(Nblock);
        J(l, c) = randn;
    end
end

F = [];
for l = 1 : numBlockL
    
    FF = [];
    for c = 1 : numBlockC
        mat = dftmtx(Nblock);
        mmat = J(l, c) * mat(rp{l, c}(1 : Mblock(l) ), :);
        mmat(noBlockError{l}, :) = -mmat(noBlockError{l}, :);
        FF = [FF, mmat];
    end
    F = [F; FF];
end

% IF = [];
% for c = 1 : numBlockC
%
%     IFF = [];
%     for l = 1 : numBlockL
%         mat = conj(dftmtx(Nblock) ) .* J(l, c);
%         mmat = zeros(Nblock);
%         mmat(:, rp{l, c}(1 : Mblock(l) ) ) = mat(:, rp{l, c}(1 : Mblock(l) ) );
%         mmat(:, noBlockError{l} ) = -mmat(:, noBlockError{l} );
%         IFF = [IFF, mmat];
%     end
%     IF = [IF; IFF];
% end

IF = [];
for c = 1 : numBlockC
    
    IFF = [];
    for l = 1 : numBlockL
        mat = conj(dftmtx(Nblock) ) .* J(l, c);
        mmat = mat(:, rp{l, c}(1 : Mblock(l) ) );
        mmat(:, noBlockError{l} ) = -mmat(:, noBlockError{l} );
        IFF = [IFF, mmat];
    end
    IF = [IF; IFF];
end

% SMtoN = zeros(1, Nblock);
% SMtoN(rp{1,1}(1 : Mblock(1) ) ) = SM;

F_SN = MultSeededFourier(SN, J, numBlockL, numBlockC, Mblock, Nblock, rp, noBlockError);
F_SN2 = MultSeededFourierSquarred(SN, J, numBlockL, numBlockC, Mblock, Nblock);

IF_SM = MultSeededFourierInverse(SM, J, numBlockL, numBlockC, Mblock, Nblock, rp, noBlockError);
IF_SM2 = MultSeededFourierInverseSquarred(SM, J, numBlockL, numBlockC, Mblock, Nblock);

max(abs(real(F * SN) - real(F_SN) ) )
max(abs(imag(F * SN) - imag(F_SN) ) )

max(abs(real(abs(F).^2 * SN) - real(F_SN2) ) )
max(abs(imag(abs(F).^2 * SN) - imag(F_SN2) ) )

max(abs(real(IF * SM) - real(IF_SM) ) )
max(abs(imag(IF * SM) - imag(IF_SM) ) )

max(abs(real(abs(IF).^2 * SM) - real(IF_SM2) ) )
max(abs(imag(abs(IF).^2 * SM) - imag(IF_SM2) ) )