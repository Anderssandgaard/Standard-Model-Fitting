function [pars,sse,flag,n] = comps_SM(data,pars,model_info,fit_info)

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
Dlb = 0.1;
Dub = 2.9;

%    [S1 Dapar S2  Depar Deperp] %Kernel parametre
lb = [0   0.01  0 0.01 0.01 ];
ub = [inf 3   inf 3 3];

lb = lb.';
ub = ub.';
fit_info.lb = lb;
fit_info.ub = ub;

N = size(pars(:,:),2);
dims = size(pars);
dims(end+1) = 1;
sse = zeros(dims(2:end),'like',data);
flag = zeros(dims(2:end),'like',data);
n = zeros(dims(2:end),'like',data);
if fit_info.display
    WaitMessage = levmar.parfor_wait(N, 'Waitbar', true);
else
    WaitMessage = [];
end
parfor i = 1:N
    y = data(:,i);
    [pars(:,i),sse(i),flag(i),n(i)] = levmar.levmar_2step_varproj(@(x)fun.comps_SM(x,y,model_info),pars(:,i),y,fit_info);
    if fit_info.display
        WaitMessage.Send;
    end
end
if fit_info.display
    levmar.warnings(flag)
end

