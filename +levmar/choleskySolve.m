function b = choleskySolve(A,b)
% solves linear system Ax=b by cholesky decomposition. A must be symmetric
% and positive-definite. The solution is in-place and goes into b. A is
% destroyed during the Cholesky decomposition.

% cholesky decomposition
decomp = levmar.cholesky(A);

% solve Ly=b by forward substitution
for m = 1:length(b)
    %for i = 1:m-1
    %    b(m) = b(m) - decomp(m,i)*b(i);
    %end
    b(m) = b(m) - decomp(m,1:m-1)*b(1:m-1);
end

% solve DL'x=y by back substitution
for m = fliplr(1:length(b))
    b(m) = b(m)/decomp(m,m);
    %for i = m+1:length(b)
    %    b(m) = b(m) - decomp(i,m)*b(i);
    %end
    b(m) = b(m) - decomp(m+1:end,m).'*b(m+1:end);
end