function [x,sse,flag,n] = fit_class(obj,y,x,imax,gradeps,display)
if nargin<5 || isempty(gradeps)
    gradeps = 0;
end
if nargin<6 || isempty(display)
    display = true;
end

N = size(x(:,:),2);
dims = size(x);
dims(end+1) = 1;
sse = zeros(dims(2:end));
flag = zeros(dims(2:end));
n = zeros(dims(2:end));
lambda = 1e-3;
upFactor = 3;
downFactor = 2;
if display
    WaitMessage = levmar.parfor_wait(N, 'Waitbar', true);
else
    WaitMessage = [];
end
parfor i = 1:N
    [x(:,i),sse(i),flag(i),n(i)] = levmar.levmar_class(obj,x(:,i),y(:,i),obj.lb,obj.ub,lambda,upFactor,downFactor,gradeps,imax);
    if display
        WaitMessage.Send;
    end
end