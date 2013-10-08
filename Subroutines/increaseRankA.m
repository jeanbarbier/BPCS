rank(A)
for i = 1 : 1000
    col = ceil(rand * numBlockC);
    AA = A;
    l1_ = ceil(N*rand);
    l2_ = ceil(N*rand);
    l1 = AA(l1_,Nblock * (col - 1) + 1 : Nblock * col);
    l2 = AA(l2_,Nblock * (col - 1) + 1 : Nblock * col);
    AA(l1_,Nblock * (col - 1) + 1 : Nblock * col) = l2;
    AA(l2_,Nblock * (col - 1) + 1 : Nblock * col) = l1;
    A = AA;
end
rank(A)
