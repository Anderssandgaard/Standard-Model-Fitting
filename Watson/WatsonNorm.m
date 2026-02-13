function p = WatsonNorm(scale,DirWatson)

b = @(theta,phi) exp((cos(phi).*sin(theta)*DirWatson(1)+...
    sin(phi).*sin(theta)*DirWatson(2)+cos(theta)*DirWatson(3)).^2*scale).*sin(theta);

p = integral2(b,0,pi,0,2*pi);
