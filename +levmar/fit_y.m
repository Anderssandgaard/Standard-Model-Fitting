function [x,sse,flag,n] = fit_y(x,y,fun,gradeps,imax,lb,ub,lambda,upFactor,downFactor,display,minRelDx)
if nargin<12
    minRelDx = 0;
end
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
for i = 1:N
    [x(:,i),sse(i),flag(i),n(i)] = levmar.levmar_y(fun,x(:,i),y(:,i),lb,ub,lambda,upFactor,downFactor,gradeps,imax,minRelDx);
    if display
        WaitMessage.Send;
    end
end