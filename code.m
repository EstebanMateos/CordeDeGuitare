n=power(2,8);
L = 1.5;
c = 16;
tf = L/c;
CFL = 0.8;
dx = L/(n-1);
dt = CFL*dx/c;
m = fix(tf/dt);
xmax = n;
tmax = m;
pert_init = fix(n/2);
M = zeros(n,m);

I=[2:n-1]; Ip1=[3:n]; Im1=[1:n-2];
X=[0:dx:L];


L1 = fix(n/4);



%Profil initial (Modifier valeur de PrI pour changer de profil initial)
PrI = 6;
if (PrI==1) 
  % Créneau
  delta=L/8;
  delta0=L/2;
  W=(X-delta0<=delta).*(X-delta0>=-delta); 
end
if (PrI==2)
  % Sinusoïdal
  f=220;
  W=cos(2*pi*f*X/L);
end

if (PrI==3)
% Dirac
M(pert_init,1) = 0.05;


M(pert_init,1) = 0.05;
M(pert_init+1,1) = 0.025;
M(pert_init-1,1) = 0.025;
M(pert_init+2,1) = 0.0125;
M(pert_init-2,1) = 0.0125;
M(pert_init+3,1) = 0.00625;
M(pert_init-3,1) = 0.00625;

end
if(PrI==4)
%Profil triangulaire isocèle

for i = 1:fix(n/2)
    M(i,1) = 0.05*i/fix(n/2);
end
for i = xmax:-1:fix(n/2)+1
    M(i,1) = M(n-i+1,1);
end
end

if(PrI==5)
%Profil triangulaire (sur un côté)
for i = 1:L1
    M(i,1) = 0.05*i/(L1-1);
end
for i = 1:xmax-L1
    M(i+L1,1) = M(L1,1)- i*0.05/(xmax-L1-1);
end
end

if(PrI==6)
%Profil sinusoïdal

for i=2:xmax-1
    M(i,1) = sin(2*pi*i*dx*2/L)
end

end

if(PrI==1 | PrI==2)
W(Im1)'
M(:,1) = W(1:xmax)';
end


dx

%Première itération
for i=2:xmax-1
    M(i,2) = M(i,1) + CFL^2/2*(M(i+1,1)-2*M(i,1)+M(i-1,1));
end
M(1,:) = 0;
M(xmax,:) = 0;

for j=3:tmax
    for i=2:xmax-1
       M(i,j) = CFL^2 * (M(i+1,j-1) - 2*M(i,j-1) + M(i-1,j-1)) - M(i,j-2) + 2 * M(i,j-1);
    end
end



x = 1:tmax;
y = 1:xmax;
[X,Y] = meshgrid(x,y);
mesh(X,Y,M)
colorbar
xlabel('Temps[dt]')
ylabel('Position[dx]')
zlabel('Perturbation[m]')
title('Propagation d une perturbation le long d une corde dans le temps(m)')
soundsc(M(:,n))
