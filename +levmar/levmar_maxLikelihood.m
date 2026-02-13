function [x,L,flag,n] = levmar_maxLikelihood(f,x,lb,ub,lambda,upFactor,downFactor,gradeps,imax)
% [L,J] = f(x) returns L "negative log likelihood", J "scaled Jacobian dPi/dxj/Pi"
flag = 0;
minDiagVals = 1e-3;
maxD = minDiagVals*ones(length(x),1);

x(x>ub) = ub(x>ub);
x(x<lb) = lb(x<lb);

x = x(:);
[L,Jp,JJp] = f(x);
for n = 1:imax
    % need stopping condition
    
    D = diag(JJp);
    filter = D>maxD;
    maxD(filter) = D(filter);
    
    dx = levmar.choleskySolve(JJp+lambda*diag(maxD),Jp);
    
    xtrial = x+dx;
    xtrial(xtrial>ub) = ub(xtrial>ub);
    xtrial(xtrial<lb) = lb(xtrial<lb);
    if all(xtrial==x) % overflow. Steps too small
        flag = 1;
        return
    end
    
    Ltrial = f(xtrial);
    if Ltrial<L
        x = xtrial;
        L = Ltrial;
        [~,Jp,JJp] = f(x);
        lambda = lambda/downFactor;
    else
        lambda = lambda*upFactor;
    end
end

if n==imax
    flag = 2;
end
end