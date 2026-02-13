function pl = watson_pl(k,lmax)
%DawF = DawsonF(sqrt(k));
 DawF = dawson(sqrt(k));
pl    = zeros(lmax/2,1);
pl(1) = 1/4*(3./DawF./sqrt(k) - 2 - 3./k);
if lmax~=2
    pl(2) = 1/32./k.^2.*(105 + 12*k.*(5+k) + 5*sqrt(k).*(2*k-21)./DawF);
    if lmax~=4
        pl(3) = 1/128./k.^3.*(-5*(8*k.^3+84*k.^2+378*k+693) + 21./DawF.*sqrt(k).*(4*k.^2-20*k+165));
        if lmax~=6
            pl(4) = 1/2048./k.^4.*(35*(16*k.^4+288*k.^3+2376*k.^2+10296*k+19305) + 3./DawF.*sqrt(k).*(248*k.^3-7700*k.^2+30030*k-225225));
        end
    end
end
end
