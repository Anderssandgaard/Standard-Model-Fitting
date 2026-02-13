function [k1, k0,...
          dkdS1,dkdS2, dkdDa,dkdR2a ,dkdR2Sa, dkdLama,...
          dkdDepar, dkdDeperp, dkdR2e, dkdR2Se, dkdLame,...
          dkdS10,dkdS20, dkdDa0, dkdR2a0, dkdR2Sa0, dkdLama0,...
          dkdDepar0, dkdDeperp0, dkdR2e0, dkdR2Se0, dkdLame0]...
           = SMPR_ISO_kernel(w,Ymat,pars,b,dotprod,TE,DTE,DTE2)

S1         = pars(1);
Da         = pars(2);
R2a        = pars(3);
R2Sa       = pars(4);
Lama       = pars(5);

S2         = pars(6);
Depar      = pars(7);
Deperp     = pars(8);
R2e        = pars(9);
R2Se       = pars(10);
Lame       = pars(11);

c1         = exp(-Da*dotprod-R2a*TE-1i*Lama*DTE-R2Sa*DTE2);
c2         = exp(-b*Deperp-(Depar-Deperp)*dotprod-R2e*TE-1i*Lame*DTE-R2Se*DTE2);

k          = S1*c1 + S2*c2;

k0         = 1/(4*pi)*(w'*k).';
k1         = (Ymat*(w.*k)).';


%derivatives
if nargout>2

dkdS1       = c1;
dkdS2       = c2;

dkdDa      = -dotprod*S1.*c1;
dkdR2a     = -TE*S1.*c1;
dkdR2Sa    = -DTE2*S1.*c1;
dkdLama    = -1i*DTE*S1.*c1;

dkdDepar   = -S2.*(dotprod).*c2;
dkdDeperp  = -S2.*(b-dotprod).*c2;
dkdR2e     = -TE*S2.*c2;
dkdR2Se    = -DTE2*S2.*c2;
dkdLame    = -1i*DTE*S2.*c2;


dkdS10     = 1/(4*pi)*(w'*dkdS1).';
dkdS20     = 1/(4*pi)*(w'*dkdS1).';

dkdDa0     = 1/(4*pi)*(w'*dkdDa).';
dkdR2a0    = 1/(4*pi)*(w'*dkdR2a).';
dkdR2Sa0   = 1/(4*pi)*(w'*dkdR2Sa).';
dkdLama0   = 1/(4*pi)*(w'*dkdLama).';

dkdDepar0  = 1/(4*pi)*(w'*dkdDepar).';
dkdDeperp0 = 1/(4*pi)*(w'*dkdDeperp).';
dkdR2e0    = 1/(4*pi)*(w'*dkdR2e).';
dkdR2Se0   = 1/(4*pi)*(w'*dkdR2Se).';
dkdLame0   = 1/(4*pi)*(w'*dkdLame).';


dkdS1      = (Ymat*(w.*dkdS1)).';
dkdS2      = (Ymat*(w.*dkdS2)).';

dkdDa      = (Ymat*(w.*dkdDa)).';
dkdR2a     = (Ymat*(w.*dkdR2a)).';
dkdR2Sa    = (Ymat*(w.*dkdR2Sa)).';
dkdLama    = (Ymat*(w.*dkdLama)).';

dkdDepar   = (Ymat*(w.*dkdDepar)).';
dkdDeperp  = (Ymat*(w.*dkdDeperp)).';
dkdR2e     = (Ymat*(w.*dkdR2e)).';
dkdR2Se    = (Ymat*(w.*dkdR2Se)).';
dkdLame    = (Ymat*(w.*dkdLame)).';


end
%%

