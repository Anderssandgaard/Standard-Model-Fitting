function [S,J,p] = comps_SMPR_ISO(pars,y,info,p)

%%%%%%%%%%%%%%%%%%
% (c) BY ANDERS DYHR SANDGAARD June 2024:
% Standard model fitting using integration on lebedev grid.
% here we include both susceptibility induced relaxation (R2 and R2*) and phase
% with orientation dependence for intra-axonal water.
%%%%%%%%%%%%%%%%%%

persistent K K0 dKdS dKdS0 dKdDa  dKdDa0  
persistent dKdR2Sa  dKdR2Sa0 dKdR2a  dKdR2a0 
persistent dKdsusc dKdsusc0
persistent info0
persistent yp pp

if nargout~=2
    yp    = y;
    info0 = info;

    [K, K0,...
     dKdS, dKdDa, dKdsusc, dKdR2Sa, dKdR2a,...
     dKdS0, dKdDa0, dKdsusc0, dKdR2Sa0, dKdR2a0]...
        = Kernels.Intra.SMPR_ISO_kernel(info0.nw,info0.Ymat',pars,info0.bgn,info0.dTE,info0.TE);

    K   = fun.complexcat(1,(K));
    K0  = fun.complexcat(1,(K0));

    if nargin~=5
         if info0.gpu
             p = (((K'*K)\(K'*(yp-K0))));
         else
             p = (linsolve((K'*K),(K'*(y-K0))));
        end
    end

    pp = p;
    S = (K0+K*pp);
    J = [];

else

    S = [];
    J          = zeros(size(info0.bgn,2),5+size(info0.Ymat,2),'like',info0.bgn);
    J(:,1)     = dKdS0+dKdS*pp;
    J(:,2)     = dKdDa0+dKdDa*pp;
    J(:,3)     = dKdsusc0+dKdsusc*pp;
    J(:,4)     = dKdR2a0+dKdR2a*pp;
    J(:,5)     = dKdR2Sa0+dKdR2Sa*pp;
    J          = fun.complexcat(1,J);
    J(:,6:end) = K;
end
