function [pars,sse,flag,n] = comps_SMPR_ISO(data,pars,model_info,fit_info)

%%%%%%%%%%%%%%%%%%
% ADDED BY ANDERS DYHR SANDGAARD:
%%%%%%%%%%%%%%%%%%

if ~exist('fit_info','var') || isempty(fit_info)
    fit_info = fit.options();
end

fit_info.gpu = isgpuarray(data);
model_info.gpu = isgpuarray(data);

Dlb = 0;
Dub = 3;

lb = [0   Dlb -inf 0    0 ];
ub = [inf Dub  inf inf  inf];
lb = lb.';
ub = ub.';
fit_info.lb = lb;
fit_info.ub = ub;

N = size(pars(:,:),2);
dims = size(pars);
dims(end+1) = 1;

sse  = zeros(dims(2:end),'like',data);
flag = zeros(dims(2:end),'like',data);
n    = zeros(dims(2:end),'like',data);

w               = waitbar(0,'Percent done ...');
Dm              = parallel.pool.DataQueue;

afterEach(Dm,@parforWaitbar);
parforWaitbar(w,N)

for i = 1:N
    y = data(:,i);
    [pars(:,i),sse(i),flag(i),n(i)] = levmar.levmar_2step_varproj(@(x)fun.Intra.comps_SMPR_ISO(x,y,model_info),pars(:,i),y,fit_info);
    send(Dm,[i]);
end
if fit_info.display
    levmar.warnings(flag)
end

delete(w)