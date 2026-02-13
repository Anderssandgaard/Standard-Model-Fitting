function [S,J,p] = comps_SMPR_ISO(pars,y,info,p)

%%%%%%%%%%%%%%%%%%
% (c) BY ANDERS DYHR SANDGAARD July 2022:
% Standard model fitting using integration on lebedev grid.
% here we include both susceptibility induced relaxation (R2 and R2*) and phase
% with orientation dependence.
%%%%%%%%%%%%%%%%%%

persistent dKdS1 dKdS10 dKdS2 dKdS20 K K0 dKdDa dKdDepar dKdDeperp dKdDa0 dKdDepar0 dKdDeperp0
persistent dKdR2Sa dKdR2Se dKdR2Sa0 dKdR2Se0
persistent dKdR2a dKdR2e dKdR2a0 dKdR2e0
persistent dKdLama dKdLame dKdLama0 dKdLame0 
persistent info0
persistent yp pp

if nargout~=2
    yp    = y;
    info0 = info;

    [K, K0,...
          dKdS1,dKdS2, dKdDa,dKdR2a ,dKdR2Sa, dKdLama,...
          dKdDepar, dKdDeperp, dKdR2e, dKdR2Se, dKdLame,...
          dKdS10,dKdS20, dKdDa0, dKdR2a0, dKdR2Sa0, dKdLama0,...
          dKdDepar0, dKdDeperp0, dKdR2e0, dKdR2Se0, dKdLame0]...
           = Kernels.SMPR_ISO_kernel(info0.nw,info0.Ymat',pars,info0.b,info0.bgn,info0.TE,info0.DTE,info0.DTE2);


    K         = fun.complexcat(1,(K));
    K0        = fun.complexcat(1,(K0));

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
    J          = zeros(length(info0.b),11+size(info0.Ymat,2),'like',info0.b);
    J(:,1)     = dKdS10+dKdS1*pp; 
    J(:,2)     = dKdDa0+dKdDa*pp;
    J(:,3)     = dKdR2a0+dKdR2a*pp;
    J(:,4)     = dKdR2Sa0+dKdR2Sa*pp;
    J(:,5)     = dKdLama0+dKdLama*pp;

    J(:,6)     = dKdS20+dKdS2*pp; 
    J(:,7)     = dKdDepar0+dKdDepar*pp;
    J(:,8)     = dKdDeperp0+dKdDeperp*pp;
    J(:,9)     = dKdR2e0+dKdR2e*pp;
    J(:,10)    = dKdR2Se0+dKdR2Se*pp;
    J(:,11)    = dKdLame0+dKdLame*pp;

    J          = fun.complexcat(1,J);
    J(:,12:end) = K;
end