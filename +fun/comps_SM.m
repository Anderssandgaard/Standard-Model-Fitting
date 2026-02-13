function [S,J,p] = comps_SM(pars,y,info,p)

%%%%%%%%%%%%%%%%%%
% BY ANDERS DYHR SANDGAARD July 2022:
% Standard model fitting using integration on lebedev grid.
%%%%%%%%%%%%%%%%%%

persistent K K0 dKdS1 dKdS2 dKdDa dKdDepar dKdDeperp dKdS10 dKdS20 dKdDa0 dKdDepar0 dKdDeperp0
persistent info0 
persistent yp pp

if nargout~=2
    yp    = y;
    info0 = info;

[K, K0,...
          dKdS1,dKdS2, dKdDa,...
          dKdDepar, dKdDeperp,...
          dKdS10,dKdS20, dKdDa0,...
          dKdDepar0, dKdDeperp0]...
          = Kernels.SM_kernel(info0.nw,info0.Ymat',pars,info0.b,info0.bgn);

    if nargin~=4
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
    J          = zeros(size(info0.b,2),5+size(info0.Ymat,2),'like',info0.b);
    J(:,1)     = dKdS10+dKdS1*pp; 
    J(:,2)     = dKdDa0+dKdDa*pp; 
    J(:,3)     = dKdS20+dKdS2*pp; 
    J(:,4)     = dKdDepar0+dKdDepar*pp; 
    J(:,5)     = dKdDeperp0+dKdDeperp*pp;
    J(:,6:end) = K;
end
