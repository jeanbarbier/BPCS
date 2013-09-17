function IM=haar_transform(size,num,IM)
    mysize=size;
    for iteration=1:num      
        MYPIC=IM(1:mysize,1:mysize);
        if (iteration>1)
            UPRIGHT=IM(1:mysize,mysize+1:size);
            LOWER=IM(mysize+1:size,:);
        end
        W=sparse_haar2(mysize);
        WT=W';
        TRANS=W*MYPIC*WT;
        if (iteration==1)
            IM=TRANS;
        else
            IM=[TRANS UPRIGHT;LOWER];
        end
        mysize=mysize/2;
    end
end