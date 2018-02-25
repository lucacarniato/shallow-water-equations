% This script solves the shallow water equations in 2d
% with an explicit time scheme
clc;
clear;
close all;

%domain
m       = 30; 
dx      = 1/m; 
dy      = 1/m;
xc       = -dx/2:dx:1+dx/2; 
yc       = -dy/2:dy:1+dy/2; 
[xxc,yyc] = meshgrid(xc,yc);
g       = 9.81; 
c       = 0.5;
nh      = size(xxc,1)-2; %with no ghost nodes
nhg     = size(xxc,1); %with ghost nodes 
dbug    = 0;
dt      = 0.001;
tfinal  = 3;

%initialize h, u, v
t      = 0;
h      = zeros(nhg);
h(:,:) = 5;
dspl = 1/10;
h(xxc>=0.5 & xxc<=0.6 & yyc>=0.5 & yyc<=0.6) = 10;
u      = zeros(nhg, nhg-1); 
v      = zeros(nhg-1, nhg); 
uup    = u;
vup    = v;
unew   = u;
vnew   = v;
hnew   = h;
uflux  = zeros(nh,nh);
vflux  = zeros(nh,nh);
%matrix diagonals
diagh  = zeros(nhg*nhg,5);
rhsh   = zeros(nhg*nhg,1);
diagu  = zeros(size(u,1) * size(u,2),5);
rhsu   = zeros(size(u,1) * size(u,2),1);
diagv  = zeros(size(v,1) * size(v,2),5);
rhsv   = zeros(size(v,1) * size(v,2),1);
%solutions
unewvec = rhsu;
uoldvec = rhsu;
vnewvec = rhsv;
voldvec = rhsv;
holdvec = rhsh;
hnewvec = rhsh;
maxflux =0;
for i=1:nhg
diaghind = nhg * (i-1) + 1 : nhg * i;
hnewvec(diaghind) = h(i,:);    
end
nf = 0;
while ( t < tfinal)
%here start the time loop
%assign new solution
h       = hnew;
u       = unew;
v       = vnew;
holdvec = hnewvec;
uoldvec = unewvec;
voldvec = vnewvec;
uup     = u >= 0; %upwind weight
vup     = v >= 0; %upwind weight
%plot solution
mesh(xxc(2:nhg,2:nhg),yyc(2:nhg,2:nhg),h(2:nhg,2:nhg)), colormap gray, axis([0 1 0 1 0.5 12]);
title(['t = ' num2str(t)])
xlabel x, ylabel y; zlabel h; pause(0.001);
nf=nf+1;
%    
for i=2:nhg-1
%compute flux, for later dt computation    
uflux(i-1,:) = 0.5 * abs(uup(i,2:nh+1).*u(i,2:nh+1).*h(i,2:nh+1)+(1 - uup(i,2:nh+1)).*u(i,2:nh+1).*h(i,3:nh+2)+ ...
                         uup(i,1:nh  ).*u(i,1:nh  ).*h(i,1:nh  )+(1 - uup(i,1:nh  )).*u(i,1:nh ).*h(i,2:nh+1 ))./h(i,2:nh+1)+...
                         sqrt(g*h(i,2:nh+1));
%compute flux, for later dt computation
vflux(:,i-1) = 0.5 * abs(vup(2:nh+1,i).*v(2:nh+1,i).*h(2:nh+1,i)+(1 - vup(2:nh+1,i)).*v(2:nh+1,i).*h(3:nh+2,i)+ ...
                         vup(1:nh,i  ).*v(1:nh,i  ).*h(1:nh,i  )+(1 - vup(1:nh,i  )).*v(1:nh,i ).*h(2:nh+1,i ))./h(2:nh+1,i)+...
                         sqrt(g*h(2:nh+1,i));
end
maxflux = max([max(max(uflux)) max(max(vflux))]);
dt      = dx *0.7/maxflux;
ntsx    = dt/dx; %normalized timestep
ntsy    = dt/dy; %normalized timestep
%--------------------------------------------------------------------------%
% continuity equation
%--------------------------------------------------------------------------%
% u component, row wise
for i=2:nhg-1
diaghind            = nhg * (i-1) + 2 : nhg * (i-1) + nhg - 1;
diagh(diaghind-1,2) = ntsx * -uup(i,1:nh).*u(i,1:nh);
diagh(diaghind,3)   = ntsx * (uup(i,2:nh+1).*u(i,2:nh+1) -(1 - uup(i,1:nh)).*u(i,1:nh));
diagh(diaghind+1,4) = ntsx * (1 - uup(i,2:nh+1)).*u(i,2:nh+1);
end
% v component, column wise
for i=2:nhg-1
diaghind              = nhg + i : nhg : nhg * (nhg-1); %strided access..
diagh(diaghind-nhg,1) = ntsy * -vup(1:nh,i).*v(1:nh,i);
diagh(diaghind,3)     = diagh(diaghind,3)   + ntsy * (vup(2:nh+1,i).*v(2:nh+1,i) -(1 - vup(1:nh,i)).*v(1:nh,i));
diagh(diaghind+nhg,5) = ntsy * (1 - vup(2:nh+1,i)).*v(2:nh+1,i);
end
hmat  = spdiags(diagh,[-nhg -1 0 1 nhg],nhg*nhg,nhg*nhg);
if (dbug)
fullhmat = full(hmat);
end
%Explicit solution
hnewvec = holdvec - hmat*holdvec;
for i=1:nhg
vecind = nhg * (i-1) + 1 : nhg * i;
hnew(i,:)= hnewvec(vecind);    
end
%--------------------------------------------------------------------------%
%U MOM
%--------------------------------------------------------------------------%
for i=2:nhg-1
diaghuind  = ((nh+1) * (i-1) + 2) : ((nh+1) * (i-1) + nh);

intup  = ntsx *(0.5 * (uup(i,3:nh+1).* h(i,3:nh+1)   + (1-uup(i,3:nh+1)).* h(i,4:nh+2)).*u(i,3:nh+1) + ...
                0.5 * (uup(i,2:nh  ).* h(i,2:nh  )   + (1-uup(i,2:nh  )).* h(i,3:nh+1)).*u(i,2:nh  )); 
    
intum  = ntsx *(0.5 * (uup(i,2:nh  ).* h(i,2:nh  )   + (1-uup(i,2:nh  )).* h(i,3:nh+1)).*u(i,2:nh  ) + ...
                0.5 * (uup(i,1:nh-1).* h(i,1:nh-1)   + (1-uup(i,1:nh-1)).* h(i,2:nh  )).*u(i,1:nh-1)); 
    
intvup = ntsy *(0.5 * (vup(i,3:nh+1).* h(i,3:nh+1) + (1-vup(i,3:nh+1)).* h(i+1,3:nh+1)).*v(i,3:nh+1) + ...
                0.5 * (vup(i,2:nh).*   h(i,2:nh)   + (1-vup(i,2:nh  )).* h(i+1,2:nh  )).*v(i,2:nh ));
     
intvum = ntsy *(0.5 * (vup(i-1,3:nh+1).* h(i-1,3:nh+1) + (1-vup(i-1,3:nh+1)).* h(i,3:nh+1)).*v(i-1,3:nh+1) + ...
                0.5 * (vup(i-1,2:nh  ).* h(i-1,2:nh  ) + (1-vup(i-1,2:nh)  ).* h(i,2:nh)  ).*v(i-1,2:nh));

havold = 0.5 * (h(i,2:nh ) + h(i,3:nh+1));
   
intupp  = intup>=0;
intump  = intum>=0;
intvupp = intvup>=0;
intvump = intvum>=0;

diagu(diaghuind-(nh+1),1) = -intvum.*(intvump); 
diagu(diaghuind-1,2)      = -intum.* intump;
diagu(diaghuind,3)        = - havold + intup.* intupp - intum.*(1 - intump) + intvup.*intvupp - intvum.*(1 - intvump);
diagu(diaghuind+1,4)      =  intup.* (1-intupp);
diagu(diaghuind+(nh+1),5) =  intvup.* (1-intvump);

%gravity rhs
rhsu(diaghuind) = ntsx * 0.5 * g * (hnew(i,3:nh+1).^2 - hnew(i,2:nh).^2);

end
% matrix assembly 
umat  = spdiags(diagu,[-(nh+1) -1 0 1 (nh+1)],size(diagu,1),size(diagu,1));
if (dbug)
fullumat = full(umat);
end 
unewvec = -umat * uoldvec - rhsu;
ind = 0;
for i=1:size(u,1)
    for j=1:size(u,2)
        ind      = ind + 1;
        havnew   = 0.5 * (hnew(i,j ) + hnew(i,j+1));
        if (havnew>1e-12)
        unewvec(ind) = unewvec(ind)/havnew;
        unew(i,j) = unewvec(ind);
        else
        unewvec(ind) =0 ;
        unew(i,j) = 0;
        end 
    end
end
%--------------------------------------------------------------------------%
%V MOM
%--------------------------------------------------------------------------%
for i=2:nhg-1
diaghvind  = ((nh+1) * (i-1) + 2 ): ((nh+1) * (i-1) + nh);

intvp  = ntsy *(0.5 * (vup(3:nh+1,i).* h(3:nh+1,i)   + (1-vup(3:nh+1,i)).* h(4:nh+2,i)).*v(3:nh+1,i) + ...
                0.5 * (vup(2:nh  ,i).* h(2:nh,i  )   + (1-vup(2:nh  ,i)).* h(3:nh+1,i)).*v(2:nh  ,i)); 
    
intvm  = ntsy *(0.5 * (vup(2:nh,i  ).* h(2:nh,i  )   + (1-vup(2:nh,i  )).* h(3:nh+1,i)).*v(2:nh,i  ) + ...
                0.5 * (vup(1:nh-1,i).* h(1:nh-1,i)   + (1-vup(1:nh-1,i)).* h(2:nh,i  )).*v(1:nh-1,i)); 
    
intuvp = ntsx *(0.5 * (uup(3:nh+1,i).* h(3:nh+1,i) + (1-uup(3:nh+1,i)).* h(3:nh+1,i+1)).*u(3:nh+1,i) + ...
                0.5 * (uup(2:nh ,i).*  h(2:nh,i)   + (1-uup(2:nh,i)).*   h(2:nh,i+1)).*  u(2:nh,i));
     
intuvm = ntsx *(0.5 * (uup(3:nh+1,i-1).* h(3:nh+1,i-1) + (1-uup(3:nh+1,i-1)).* h(3:nh+1,i)).*u(3:nh+1,i-1) + ...
                0.5 * (uup(2:nh,i-1).*   h(2:nh,i-1  ) + (1-uup(2:nh,i-1  )).* h(2:nh,i  )).*u(2:nh,i-1  ));

havold = 0.5 * (h(2:nh,i) + h(3:nh+1,i));
  
intvpp  = intvp>=0;
intvmp  = intvm>=0;
intuvpp = intuvp>=0;
intuvmp = intuvm>=0;

diagv(diaghvind-(nh+1),1) = -intuvm.*(intuvmp); 
diagv(diaghvind-1,2)      = -intvm.* intvmp;
diagv(diaghvind,3)        = -havold + intvp.* intvpp - intvm.*(1 - intvmp) + intuvp.*intuvpp - intuvm.*(1 - intuvmp);
diagv(diaghvind+1,4)      =  intvp.* (1-intvpp);
diagv(diaghvind+(nh+1),5) =  intuvp.* (1-intuvmp);

%gravity rhs
rhsv(diaghvind) = ntsy * 0.5 * g * (hnew(3:nh+1,i).^2 - hnew(2:nh,i).^2);
end
% matrix assembly  
vmat  = spdiags(diagv,[-(nh+1) -1 0 1 (nh+1)],size(diagv,1),size(diagv,1));
if (dbug)
fullvmat = full(vmat);
end
vnewvec = -vmat * voldvec - rhsv;
for i=1:size(v,1)
    for j=1:size(v,2)
        ind      = i + (j-1) * size(v,1);
        havnew   = 0.5 * (hnew(i,j) + hnew(i+1,j));
        if (havnew>1e-12)
        vnewvec(ind)=vnewvec(ind)/havnew;
        vnew(i,j) = vnewvec(ind);
        else
        vnewvec(ind) =0 ;
        vnew(i,j) = 0;
        end 
    end
end
%update timestep
t = t + dt;
end