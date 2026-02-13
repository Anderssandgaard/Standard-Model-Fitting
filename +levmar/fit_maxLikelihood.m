function [x,L,flag,n] = fit_maxLikelihood(x,y,fun,gradeps,imax,lb,ub,lambda,upFactor,downFactor,display)
N = size(x(:,:),2);
dims = size(x);
dims(end+1) = 1;
L = zeros(dims(2:end));
flag = zeros(dims(2:end));
n = zeros(dims(2:end));
if display
    WaitMessage = levmar.parfor_wait(N, 'Waitbar', true);
else
    WaitMessage = [];
end
parfor i = 1:N
    [x(:,i),L(i),flag(i),n(i)] = levmar.levmar_maxLikelihood(@(x)fun(x,y(:,i)),x(:,i),lb,ub,lambda,upFactor,downFactor,gradeps,imax);
    [x(:,i),L(i)] = fminsearch(@(x)fun(x,y(:,i)),x(:,i)); % temporary
    if display
        WaitMessage.Send;
    end
end