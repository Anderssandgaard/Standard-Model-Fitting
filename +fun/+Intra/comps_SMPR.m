function [S,J,p] = comps_SMPR(pars,y,info,p)

%%%%%%%%%%%%%%%%%%
% (c) BY ANDERS DYHR SANDGAARD June 2024:
% Standard model fitting using integration on lebedev quadrature grid.
% here we include both susceptibility induced relaxation (R2 and R2*) and phase
% with orientation dependence for intra-axonal water.
%%%%%%%%%%%%%%%%%%

persistent K K0 dKdS dKdS0 dKdDa  dKdDa0  
persistent dKdR2a  dKdR2a0 
persistent dKdR2sa  dKdR2sa0 
persistent dKdLam dKdLam_ani dKdLam0 dKdLam_ani0
persistent info0 yp
persistent pp

if nargout~=2
    yp    = y;
    info0 = info;

    [K, K0,...
     dKdS, dKdDa,dKdLam,dKdLam_ani,dKdR2a,dKdR2sa,...
     dKdS0, dKdDa0,dKdLam0,dKdLam_ani0, dKdR2a0,dKdR2sa0]...
        = Kernels.Intra.SMPR_kernel(info0.nw,info0.Ymat',pars,...
                                    info0.bgn,info0.TE,...
                                    info0.dTE,info0.dTER);

    K   = fun.complexcat(1,(K));
    K0  = fun.complexcat(1,(K0));

    if nargin~=5
        if info0.gpu
            p = (((K'*K)\(K'*(yp-K0))));
        else
            p = (linsolve((K'*K),(K'*(yp-K0))));
        end
    end

    pp = p;
    S = (K0+K*pp);
    J = [];

else

    S = [];
    J           = zeros(size(info0.bgn,2),...
                   6+size(info0.Ymat,2),'like',info0.bgn);
    J(:,1)      = dKdS0+dKdS*pp;
    J(:,2)      = dKdDa0+dKdDa*pp;
    J(:,3)      = dKdLam0+dKdLam*pp;
    J(:,4)      = dKdLam_ani0+dKdLam_ani*pp;
    J(:,5)      = dKdR2a0+dKdR2a*pp;
    J(:,6)      = dKdR2sa0+dKdR2sa*pp;
    J           = fun.complexcat(1,J);
    J(:,7:end)  = K;

end
