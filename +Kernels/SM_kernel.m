function [k1, k0, dkdS1, dkdS2, dkdDa, dkdDepar, dkdDeperp, dkdS10,dkdS20, dkdDa0, dkdDepar0, dkdDeperp0]...
          = SM_kernel(w,Ymat,pars,b,dotprod)

S1      = pars(1);
Da      = pars(2);
S2      = pars(3);
Depar   = pars(4);
Deperp  = pars(5);

c1      = exp(-Da*dotprod);
c2      = exp(-b*Deperp-(Depar-Deperp)*dotprod);
k       = f.*c1+(1-f)*c2;

k0      = 1/(4*pi)*(w'*k)';
k1      = (Ymat*(w.*k))';

if nargout>2

dkdS1      = c1;
dkdS2      = c2;

dkdDa      = -dotprod*S1.*c1;

dkdDepar   = -S2.*(dotprod).*c2;
dkdDeperp  = -S2.*(b-dotprod).*c2;

dkdS10     = 1/(4*pi)*(w'*dkdS1)';
dkdS20     = 1/(4*pi)*(w'*dkdS2)';

dkdDa0     = 1/(4*pi)*(w'*dkdDa)';
dkdDepar0  = 1/(4*pi)*(w'*dkdDepar)';
dkdDeperp0 = 1/(4*pi)*(w'*dkdDeperp)';

dkdS1      = (Ymat*(w.*dkdS1))';
dkdS2      = (Ymat*(w.*dkdS2))';

dkdDa      = (Ymat*(w.*dkdDa))';

dkdDepar   = (Ymat*(w.*dkdDepar))';
dkdDeperp  = (Ymat*(w.*dkdDeperp))';
end
%%

