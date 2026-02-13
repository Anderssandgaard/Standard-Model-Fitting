function [S,J,p] = comps_SMPR(pars,y,info,p)

%%%%%%%%%%%%%%%%%%
% (c) BY ANDERS DYHR SANDGAARD June 2024:
% Standard model fitting using integration on lebedev grid.
% here we include both susceptibility induced relaxation (R2 and R2*) and phase
% with orientation dependence for intra-axonal water and extra-axonal water.
%%%%%%%%%%%%%%%%%%

persistent          K         K0 
persistent          dKdSa     dKdSe
persistent          dKdDa     dKdDepar   dKdDeperp
persistent          dKdLama_aniso  dKdLama_iso    dKdLame_aniso   dKdLame_iso
persistent          dKdR2Sa_aniso  dKdR2Sa_iso    dKdR2Se_aniso    dKdR2Se_iso
persistent          dKdR2a_iso     dKdR2a_aniso   dKdR2e_iso     dKdR2e_aniso
persistent          dKdSa0    dKdSe0
persistent          dKdDa0    dKdDepar0  dKdDeperp0
persistent          dKdLama_aniso0 dKdLama_iso0   dKdLame_aniso0  dKdLame_iso0
persistent          dKdR2a_iso0  dKdR2a_aniso0    dKdR2e_iso0     dKdR2e_aniso0
persistent          dKdR2Sa_iso0  dKdR2Sa_aniso0  dKdR2Se_iso0    dKdR2Se_aniso0
persistent          info0
persistent          yp        pp

if nargout~=2
    yp    = y;
    info0 = info;

  

[K, K0,...
          dKdSa,   dKdSe,...
          dKdDa,  dKdDepar, dKdDeperp,  ...
          dKdLama_aniso, dKdLama_iso, dKdLame_aniso, dKdLame_iso,  ...
          dKdR2Sa_aniso,  dKdR2Sa_iso,  dKdR2Se_aniso, dKdR2Se_iso, ...
          dKdR2a_iso,dKdR2a_aniso, dKdR2e_aniso, dKdR2e_iso,...
          dKdSa0, dKdSe0, ...
          dKdDa0,dKdDepar0,dKdDeperp0, ...
          dKdLama_aniso0,dKdLama_iso0,dKdLame_aniso0, dKdLame_iso0, ...
          dKdR2Sa_aniso0, dKdR2Sa_iso0, dKdR2Se_aniso0, dKdR2Se_iso0, ...
          dKdR2a_iso0,dKdR2a_aniso0, dKdR2e_aniso0, dKdR2e_iso0]...
          = Kernels.SMPR_kernel(info0.nw,info0.Ymat',pars,info0.b,info0.dotprod,info0.DTE,info0.DTE2,info0.DTES2,info0.DTE2S2,info0.DTE2S4,info0.TE,info0.TES,info0.TES4);


    K   = fun.complexcat(1,(K));
    K0  = fun.complexcat(1,(K0));

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
    %jacobian
    S = [];
    J          = zeros(length(info0.TE),17+size(info0.Ymat,2),'like',info0.TE);
    %intra
    J(:,1)     = dKdSa0+dKdSa*pp;
    J(:,2)     = dKdDa0+dKdDa*pp;
    J(:,3)     = dKdR2a_iso0+dKdR2a_iso*pp;
    J(:,4)     = dKdR2a_aniso0+dKdR2a_aniso*pp;
    J(:,5)     = dKdR2Sa_iso0+dKdR2Sa_iso*pp;
    J(:,6)     = dKdR2Sa_aniso0+dKdR2Sa_aniso*pp;
    J(:,7)     = dKdLama_iso0+dKdLama_iso*pp;
    J(:,8)     = dKdLama_aniso0+dKdLama_aniso*pp;

    %extra
    J(:,9)      = dKdSe0+dKdSe*pp;
    J(:,10)     = dKdDepar0+dKdDepar*pp;
    J(:,11)     = dKdDeperp0+dKdDeperp*pp;
    J(:,12)     = dKdR2e_iso0+dKdR2e_iso*pp;
    J(:,13)     = dKdR2e_aniso0+dKdR2e_aniso*pp;
    J(:,14)     = dKdR2Se_iso0+dKdR2Se_iso*pp;
    J(:,15)     = dKdR2Se_aniso0+dKdR2Se_aniso*pp;
    J(:,16)     = dKdLame_iso0+dKdLame_iso*pp;
    J(:,17)     = dKdLame_aniso0+dKdLame_aniso*pp;
    %fodf
    J          = fun.complexcat(1,J);
    J(:,18:end)= K;
end

