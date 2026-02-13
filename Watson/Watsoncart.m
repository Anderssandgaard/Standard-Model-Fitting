function p = Watsoncart(x1,x2,x3,scale,DirWatson,normz)

b = (x1*DirWatson(1)+x2*DirWatson(2)+x3*DirWatson(3)).^2*scale;
p = 1/normz*exp(b);
