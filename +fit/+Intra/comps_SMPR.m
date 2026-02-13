function [pars,sse,flag,n] = comps_SMPR(data,pars,model_info,fit_info)

%%%%%%%%%%%%%%%%%%
% ADDED BY ANDERS DYHR SANDGAARD:
%%%%%%%%%%%%%%%%%%

if ~exist('fit_info','var') || isempty(fit_info)
    fit_info = fit.options();
end

fit_info.gpu = isgpuarray(data);
model_info.gpu = isgpuarray(data);

lb = [];
ub = [];

lb = [0   0 -inf -inf   0     0  ];
ub = [inf 4  inf  inf   inf   inf];
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
Dm               = parallel.pool.DataQueue;

afterEach(Dm,@parforWaitbar);
parforWaitbar(w,N)

for i = 1:N
    y = data(:,i);
    send(Dm,[i]);

end
if fit_info.display
    levmar.warnings(flag)
end

delete(w)