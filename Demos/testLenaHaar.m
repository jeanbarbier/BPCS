S = double(Lena((ll - 1) .* sqrt(N) + 1: ll .* sqrt(N), (cc - 1) .* sqrt(N) + 1 : cc .* sqrt(N) ) ) ./ 256; % patch of Lena
Sh = haar_transform(sqrt(N), 5, S);
[SsortBig, IndSsortBig] = sort(reshape(Sh,1,N),'descend');
[SsortSmall, IndSsortSmall] = sort(reshape(Sh,1,N),'ascend');

Sbig = SsortBig(1 : ceil(N * 0.1) );
Ssmall = SsortSmall(1 : ceil(N * 0.1) );
Stry = zeros(1,N);
Stry(IndSsortBig(1 : ceil(N * 0.1) ) ) = Sbig;
Stry(IndSsortSmall(1 : ceil(N * 0.1) ) ) = Ssmall;

imshow(haar_inverse(sqrt(N), 5, reshape(Stry,sqrt(N),sqrt(N))))
