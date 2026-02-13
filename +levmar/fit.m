function [x,sse,flag,n] = fit(x,y,w,fun,gradeps,imax,lb,ub,lambda,upFactor,downFactor,display)
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
    [x(:,i),sse(i),flag(i),n(i)] = levmar.levmar(fun,x(:,i),y(:,i),w,lb,ub,lambda,upFactor,downFactor,gradeps,imax);
    if display
        WaitMessage.Send;
    end
end