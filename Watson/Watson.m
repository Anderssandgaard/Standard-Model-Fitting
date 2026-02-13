function p = Watson(theta,phi,scale,DirWatson,normz)

x1 = cos(phi).*sin(theta);
x2 = sin(phi).*sin(theta);
x3 = cos(theta);

b = (x1*DirWatson(1)+x2*DirWatson(2)+x3*DirWatson(3)).^2*scale;
p = 1/normz*exp(b);
