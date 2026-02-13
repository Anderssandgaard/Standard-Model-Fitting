
function [K1, K0,...
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
          = SMPR_kernel(w,Ymat,pars,b,dotprod,DTE,DTE2,DTES2,DTE2S2,DTE2S4,TE,TES,TES4)

%intra
Sa         = pars(1);
Da         = pars(2);
R2a_iso    = pars(3);
R2a_aniso  = pars(4);
R2Sa_iso   = pars(5);
R2Sa_aniso = pars(6);
Lama_iso   = pars(7);
Lama_aniso = pars(8);

%extra
Se         = pars(9);
Depar      = pars(10);
Deperp     = pars(11);
R2e_iso    = pars(12);
R2e_aniso  = pars(13);
R2Se_iso   = pars(14);
R2Se_aniso = pars(15);
Lame_iso   = pars(16);
Lame_aniso = pars(17);


% Kernel matrix
ca         = exp(-Da*dotprod...
                    -1i*Lama_aniso*DTES2-1i*Lama_iso*DTE...
                    -R2Sa_aniso*DTE2S2-R2Sa_iso*DTE2...
                    -R2a_aniso*TES-R2a_iso*TE); %intra

ce         = exp(-b*Deperp-(Depar-Deperp)*dotprod...
                    -1i*Lame_aniso*DTES2-1i*Lame_iso*DTE...
                    -R2Se_aniso*DTE2S4-R2Se_iso*DTE2...
                    -R2e_aniso*TES4-R2e_iso*TE); %extra

k          = Sa*ca + Se*ce;
K0         = 1/(4*pi)*(w'*k).';
K1         = (Ymat*(w.*k)).';

% Derivatives for jacobian
if nargout>2

%intra    
dKdSa            = ca;
dKdDa            = -dotprod*Sa.*ca;
dKdR2a_iso       = -TE*Sa.*ca;
dKdR2a_aniso     = -TES*Sa.*ca;
dKdR2Sa_iso      = -DTE2*Sa.*ca;
dKdR2Sa_aniso    = -DTE2S2*Sa.*ca;
dKdLama_iso      = -1i*DTE*Sa.*ca;
dKdLama_aniso    = -1i*DTES2*Sa.*ca;

dKdSa0           = 1/(4*pi)*(w'*dKdSa).';
dKdDa0           = 1/(4*pi)*(w'*dKdDa).';
dKdLama_iso0     = 1/(4*pi)*(w'*dKdLama_iso).';
dKdLama_aniso0   = 1/(4*pi)*(w'*dKdLama_aniso).';
dKdR2Sa_iso0     = 1/(4*pi)*(w'*dKdR2Sa_iso).';
dKdR2Sa_aniso0   = 1/(4*pi)*(w'*dKdR2Sa_aniso).';
dKdR2a_iso0      = 1/(4*pi)*(w'*dKdR2a_iso).';
dKdR2a_aniso0    = 1/(4*pi)*(w'*dKdR2a_aniso).';

dKdSa            = (Ymat*(w.*dKdSa)).';
dKdDa            = (Ymat*(w.*dKdDa)).';
dKdLama_iso      = (Ymat*(w.*dKdLama_iso)).';
dKdLama_aniso    = (Ymat*(w.*dKdLama_aniso)).';
dKdR2Sa_iso      = (Ymat*(w.*dKdR2Sa_iso)).';
dKdR2Sa_aniso    = (Ymat*(w.*dKdR2Sa_aniso)).';
dKdR2a_iso       = (Ymat*(w.*dKdR2a_iso)).';
dKdR2a_aniso     = (Ymat*(w.*dKdR2a_aniso)).';

%extra    
dKdSe            = ce;
dKdDepar         = -Se.*(dotprod).*ce;
dKdDeperp        = -Se.*(b-dotprod).*ce;
dKdR2e_iso       = -TE.*ce;
dKdR2e_aniso     = -TES4*Se.*ce;
dKdR2Se_iso      = -DTE2*Se.*ce;
dKdR2Se_aniso    = -DTE2S4*Se.*ce;
dKdLame_iso      = -1i*DTE*Se.*ce;
dKdLame_aniso    = -1i*DTES2*Se.*ce;


dKdSe0           = 1/(4*pi)*(w'*dKdSe).';
dKdDepar0        = 1/(4*pi)*(w'*dKdDepar).';
dKdDeperp0       = 1/(4*pi)*(w'*dKdDeperp).';
dKdLame_iso0     = 1/(4*pi)*(w'*dKdLame_iso).';
dKdLame_aniso0   = 1/(4*pi)*(w'*dKdLame_aniso).';
dKdR2Se_iso0     = 1/(4*pi)*(w'*dKdR2Se_iso).';
dKdR2Se_aniso0   = 1/(4*pi)*(w'*dKdR2Se_aniso).';
dKdR2e_iso0      = 1/(4*pi)*(w'*dKdR2e_iso).';
dKdR2e_aniso0    = 1/(4*pi)*(w'*dKdR2e_aniso).';

dKdSe            = (Ymat*(w.*dKdSe)).';
dKdDepar         = (Ymat*(w.*dKdDepar)).';
dKdDeperp        = (Ymat*(w.*dKdDeperp)).';
dKdR2e_iso       = (Ymat*(w.*dKdR2e_iso)).';
dKdR2e_aniso     = (Ymat*(w.*dKdR2e_aniso)).';
dKdR2Se_iso      = (Ymat*(w.*dKdR2Se_iso)).';
dKdR2Se_aniso    = (Ymat*(w.*dKdR2Se_aniso)).';
dKdLame_iso      = (Ymat*(w.*dKdLame_iso)).';
dKdLame_aniso    = (Ymat*(w.*dKdLame_aniso)).';

end
%%

