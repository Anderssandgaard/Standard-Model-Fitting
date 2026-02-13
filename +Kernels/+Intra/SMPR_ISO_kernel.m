function [k1, k0,...
          dkdS, dkdDa,dkdLam,dkdR2a,dkdR2sa,...
          dkdS0, dkdDa0,dkdLam0,dkdR2a0,dkdR2sa0]...
           = SMPR_ISO_kernel(w,Ymat,pars,dotprod,TE,DTE)


S0          = pars(1);
Da          = pars(2);
Lam         = pars(3);
R2a         = pars(4);
R2sa        = pars(5);

k           = S0*exp(-Da*dotprod... %diffusion
                     -R2a*TE... %TE part
                     -1i*Lam*DTE-R2sa*DTE); %DTE part

k0          = 1/(4*pi)*(w'*k).';
k1          = 1/(4*pi)*(Ymat*(w.*k)).';

%derivatives
if nargout>2
    
dkdS         = k/S0;
dkdDa        = -dotprod.*k;
dkdLam       = -1i*DTE.*k;
dkdR2a       = -TE.*k;
dkdR2sa      = -DTE.*k;

dkdS0        = 1/(4*pi)*(w'*dkdS).';
dkdDa0       = 1/(4*pi)*(w'*dkdDa).';
dkdLam0      = 1/(4*pi)*(w'*dkdLam).';
dkdR2a0      = 1/(4*pi)*(w'*dkdR2a).';
dkdR2sa0     = 1/(4*pi)*(w'*dkdR2sa).';

dkdS         = 1/(4*pi)*(Ymat*(w.*dkdS)).';
dkdDa        = 1/(4*pi)*(Ymat*(w.*dkdDa)).';
dkdLam       = 1/(4*pi)*(Ymat*(w.*dkdLam)).';
dkdR2a       = 1/(4*pi)*(Ymat*(w.*dkdR2a)).';
dkdR2sa      = 1/(4*pi)*(Ymat*(w.*dkdR2sa)).';

end
%%

