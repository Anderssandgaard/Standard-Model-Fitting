function [x,s2,flag,n] = levmar_varprojComplete(fun,x,Ncomps,y,lb,ub,lambda,upFactor,downFactor,gradeps,imax)
f = x(1:Ncomps);
x = x(Ncomps+1:end);
x(x>ub) = ub(x>ub);
x(x<lb) = lb(x<lb);

[S,f,J] = fun(x,f,y);
r = y - S;
s2 = r'*r;
maxD = 1e-3*ones(length(x),1);
for n = 1:imax
    J(:,x==ub) = (J(:,x==ub)<=0).*J(:,x==ub);
    J(:,x==lb) = (J(:,x==lb)>=0).*J(:,x==lb);
    Jtr = J'*r;
    if gradeps~=0 && norm(Jtr)/length(r)/length(x)<gradeps
        flag = 0;
        break
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
        break
    end
    
    [S,ftrial,Jtrial] = fun(xtrial,f,y);
    rtrial = y-S;
    s2trial = rtrial.'*rtrial;
%     [S,ftrial] = fun(xtrial,f,y);
%     s2trial = (y-S)'*(y-S);
    if s2trial<s2
        x = xtrial;
        f = ftrial;
        r = rtrial;
        s2 = s2trial;
        J = Jtrial;
%         [S,f,J] = fun(x,f,y);
%         r = y - S;
%         s2 = r'*r;
        lambda = lambda/downFactor;
    else
        lambda = lambda*upFactor;
    end
end
x = [f; x];

if n==imax
    flag = 2;
end
end