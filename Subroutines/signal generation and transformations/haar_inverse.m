function IM=haar_inverse(size,num,IM)
    mysize=size*2;
    for iteration=1:num
        mysize=mysize/2;
    end
    %Now I have the size I used!
    for iteration=1:num      
        MYPIC=IM(1:mysize,1:mysize);
        if (mysize<size)
            UPRIGHT=IM(1:mysize,mysize+1:size);
            LOWER=IM(mysize+1:size,:);
        end
        H=sparse_haar2(mysize);
        W=inv(H);
        WT=inv(H');
        TRANS=W*MYPIC*WT;
        if (mysize==size)
            IM=TRANS;
        else
            IM=[TRANS UPRIGHT;LOWER];
        end
        mysize=mysize*2;
    end
end