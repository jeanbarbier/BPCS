function [Z] = MultSeededFourierInverseSquarred(X, J, numBlockL, numBlockC, Mblock, Nblock)
% Nbloc, size of blocks, must be a power of 2
% rp*, random permutation vector for the blocks

if (isrow(X) ); X = X.'; X = [X; zeros(2 * Nblock, 1)];
else X = [X; zeros(2 * Nblock, 1)]; end
Z = zeros(numBlockC * Nblock, 1);

lastZ = 1;
for c = 1 : numBlockC
    
    u = 1;
    Y = 0;
    for l = 1 : numBlockL
        
        if (J(l, c) ~= 0); Y = Y + J(l, c)^2 * sum(X(u : u + Mblock(l) - 1) ) * ones(Nblock, 1); end;
        
        u = u + Mblock(l);
    end
    
    Z(lastZ : lastZ + Nblock - 1) = Y;
    lastZ = lastZ + Nblock;
end

if (isrow(X) ); Z = Z.'; end

end

% function [Z] = MultSeededFourierInverseSquarred(X, J, numBlockL, numBlockC, Mblock, Nblock, rp, noBlockError)
% % Nbloc, size of blocks, must be a power of 2
% % rp*, random permutation vector for the blocks
% 
% if (isrow(X) ); X = X.'; end
% Z = zeros(numBlockC * Nblock, 1);
% % ZZ = Z;
% 
% lastZ = 0; zero = zeros(1, Nblock);
% for l = 1 : numBlockL
%     
%     Y = 0; stop = Mblock(1); start = 1;
%     %     YY = 0;
%     for c = 1 : numBlockC
%         
%         if (J(l, c) ~= 0)
%             XXX = zero;
%             XXX(rp{l, c}(1 : Mblock(l) ) ) = X(start : stop);
%             XXXX = XXX(XXX ~= 0);
%             if (max(size(noBlockError) ) > 0); XXXX(noBlockError{l} ) = -XXXX(noBlockError{l} ); end
%             XXX(XXX ~= 0) = XXXX;
%             Y = Y + sum(XXX) .* J(l, c).^2;
%             %             % test
%             %                         YY = YY + (conj(dftmtx(Nblock) ) ./ Nblock .* J(l, c) ) * XXX.';
%         end
%         
%         if (l < numBlockL); start = stop + 1; stop = stop + Mblock(l); else start = stop + 1; stop = stop + Mblock(l); end
%         
%     end
%     
%     Z(lastZ + 1 : lastZ + Nblock) = Y;
%     %     ZZ(lastZ + 1 : lastZ + Nblock) = YY;
%     lastZ = lastZ + Nblock;
% end
% 
% if (isrow(X) ); Z = Z.'; end
% % max(abs(Z - ZZ) )
% 
% end