function [x,s2,flag,n] = levmar_class(obj,x,y,lb,ub,lambda,upFactor,downFactor,gradeps,imax)
x(x>ub) = ub(x>ub);
x(x<lb) = lb(x<lb);
[fx,J] = eval(obj,x,y);
r = y - fx;
s2 = r'*r;
maxD = 1e-3*ones(length(x),1);
for n = 1:imax
    J(:,x==ub) = (J(:,x==ub)<=0).*J(:,x==ub);
    J(:,x==lb) = (J(:,x==lb)>=0).*J(:,x==lb);
    Jtr = J'*r;
    if gradeps~=0 && norm(Jtr)/length(r)/length(x)<gradeps
        flag = 0;
        return
    end
    
    JtJ = J'*J;
    D = diag(JtJ);
    filter = D>maxD;
    maxD(filter) = D(filter);
    
    dx = levmar.choleskySolve(JtJ+lambda*diag(maxD),Jtr);
    
    xtrial = x+dx;
    xtrial(xtrial>ub) = ub(xtrial>ub);
    xtrial(xtrial<lb) = lb(xtrial<lb);
    if all(xtrial==x) % overflow. Steps too small
        flag = 1;
        return
    end
    
    [fx,Jtrial] = eval(obj,xtrial,y);
    rtrial = y - fx;
    s2trial = rtrial.'*rtrial;
    if s2trial<s2
        x = xtrial;
        r = rtrial;
        s2 = s2trial;
        J = Jtrial;
        lambda = lambda/downFactor;
    else
        lambda = lambda*upFactor;
    end
end
flag = 2;
end