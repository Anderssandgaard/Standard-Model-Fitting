function [x,sse,flag,n] = fit_varprojComplete(x,Ncomps,y,fun,gradeps,imax,lb,ub,lambda,upFactor,downFactor,display)
N = size(x(:,:),2);
dims = size(x);
dims(end+1) = 1;
sse = zeros(dims(2:end));
flag = zeros(dims(2:end));
n = zeros(dims(2:end));
if display
    WaitMessage = levmar.parfor_wait(N, 'Waitbar', true);
else
    WaitMessage = [];
end
parfor i = 1:N
    [x(:,i),sse(i),flag(i),n(i)] = levmar.levmar_varprojComplete(fun,x(:,i),Ncomps,y(:,i),lb,ub,lambda,upFactor,downFactor,gradeps,imax);
    if display
        WaitMessage.Send;
    end
end