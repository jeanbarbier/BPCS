function h=sparse_flo2(N)

%N=2^n; 
n1=2;
UP_init=[0 0];
DOWN_init=[-1/sqrt(2) 1/sqrt(2)];
%UP_init=[1/2 1/2];
%DOWN_init=[-1/2 1/2];
UP_matrix=[UP_init];
DOWN_matrix=[DOWN_init];
while (n1<N)
        [a,b]=size(UP_matrix);
        zero=sparse(a,b);
        UP_matrix=[UP_matrix zero;zero UP_matrix];
        DOWN_matrix=[DOWN_matrix zero;zero DOWN_matrix];
      n1=n1*2;
end
h=[UP_matrix;DOWN_matrix];