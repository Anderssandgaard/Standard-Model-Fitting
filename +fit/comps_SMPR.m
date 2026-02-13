function [pars,sse,flag,n] = comps_SMPR(data,pars,model_info,fit_info)

%%%%%%%%%%%%%%%%%%
% ADDED BY ANDERS DYHR SANDGAARD:
%%%%%%%%%%%%%%%%%%

if ~exist('fit_info','var') || isempty(fit_info)
    fit_info = fit.options();
end

fit_info.gpu   = isgpuarray(data);
model_info.gpu = isgpuarray(data);

lb    = [];
ub    = [];
Slb   = 0;
Sub   = inf;
Dlb   = 0;
Dub   = 3;
Llb   = inf;
Lub   = -inf;
R2Slb = 0;
R2Sub = inf;
R2lb  = 0;
R2ub  = inf;


lb = [Slb Dlb Llb Llb R2Slb R2Slb R2lb ... %intra
      Slb Dlb Dlb Llb Llb R2Slb R2Slb R2lb R2lb]; %extra
ub = [Sub Dub Lub Lub R2Sub R2Sub R2ub ... %intra
      Sub Dub Dub Lub Lub R2Sub R2Sub R2ub R2ub]; %extra

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

if fit_info.display
    WaitMessage = levmar.parfor_wait(N, 'Waitbar', true);
else
    WaitMessage = [];
end

for i = 1:N
    y = data(:,i);
    [pars(:,i),sse(i),flag(i),n(i)] = levmar.levmar_2step_varproj(@(x)fun.comps_SMPR(x,y,model_info),pars(:,i),y,fit_info);
    if fit_info.display
        WaitMessage.Send;
    end
end
if fit_info.display
    levmar.warnings(flag)
end

