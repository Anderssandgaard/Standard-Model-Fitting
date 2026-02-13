Angles      = [20 45 90];
Directions  = [7, 14, 21, 28];
%%
[V]=ParticleSampleSphere('N',220,'asym',true);
I = find(180/pi*acos(V(:,3))<Angles(1));
testbo = numel(I)
%%
Point_Sph.Points{1,2} = V(I,:);
Point_Sph.numsph{1,2} = 220;

%%
Point_Sph.Angles = Angles;
Point_Sph.Directions = Directions;







