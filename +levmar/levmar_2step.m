% function [x,s2,flag,n,np] = levmar_2step(f,x,y,w,lb,ub,lambda,upFactor,downFactor,gradeps,imax)
function [x,s2,flag,n,np] = levmar_2step(f,x,y,info)

lb = info.lb;
ub = info.ub;
w = info.w;
if isempty(info.w)
    w = ones(size(y));
end
imax = info.imax;
gradeps = info.gradeps;
lambda = info.lambda;
upFactor = info.upFactor;
downFactor = info.downFactor;

x(x>ub) = ub(x>ub);
x(x<lb) = lb(x<lb);
while any(x==0) % want x to be non-zero
    x(x==0) = 1e-12*randn();
    x(x>ub) = ub(x>ub);
    x(x<lb) = lb(x<lb);
end

x = x(:);
fx = f(x);
[~,J] = f(x);
r = y - fx;
s2 = r'*(r.*w);
np = 0;
maxDiag = 1e-3*ones(length(x),1);
for n = 1:imax
    J(:,x==ub) = (J(:,x==ub)<=0).*J(:,x==ub);
    J(:,x==lb) = (J(:,x==lb)>=0).*J(:,x==lb);
    Jtr = J'*(r.*w);
    if gradeps~=0 && norm(Jtr)/length(r)/length(x)<gradeps
        flag = 0;
        return
    end
    
    JtJ = J'*(J.*w);
    D = diag(JtJ);
    maxDiag(D>maxDiag) = D(D>maxDiag);
    
    dx = levmar.choleskySolve(JtJ+lambda*diag(maxDiag),Jtr);
    
    xtrial = x+dx;
    xtrial(xtrial>ub) = ub(xtrial>ub);
    xtrial(xtrial<lb) = lb(xtrial<lb);
    if all(xtrial==x) % overflow. Steps too small
        flag = 1;
        return
    end
    
    fx = f(xtrial);
    rtrial = y - fx;
    s2trial = rtrial.'*(rtrial.*w);
    if s2trial<s2
        np = np + 1;
        x = xtrial;
        r = rtrial;
        s2 = s2trial;
        [~,J] = f(x);
        if lambda>1e-6*min(maxDiag)
            lambda = lambda/downFactor;
        end
    else
        lambda = lambda*upFactor;
    end
end

flag = 2;
end