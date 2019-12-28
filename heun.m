function u=heun(n,cfl,cw,cg,T)
%int:n,gc,i
%double:cfl,cw,cg,T,g,sw_up,sw_down,sg_up,sg_down
%vector:rhoa,imua,ga,ca,b
%matrix:u0,A,u

global rhoa imua ga ca gc g;

g=1;
gc=3;

%Valores iniciales
sw_up=1.00;
sw_down=0.20;
sg_up=0.0;
sg_down=0.20;

%Constantes iniciales
imua=zeros(3,1);
imua(1)=1.0; % 1/mu
imua(2)=1.0;
imua(3)=1.0;
rhoa=zeros(3,1);
rhoa(1)=1.0;     %w
rhoa(2)=0.0012;  %g
rhoa(3)=0.85;    %o
ca=zeros(2,1);
ca(1)=cw;
ca(2)=cg;
ga=zeros(3,1);
ga(1)=0.5;
ga(2)=0.5;
ga(3)=0.5;

u0=zeros(2, n);

%Valores iniciales sw

for i=1:n/2
    u0(1, i)=sw_up;
    u0(2, i)=sg_up;
end
for i=n/2+1:n
    u0(1, i)=sw_down;
    u0(2, i)=sg_down;
end

A=zeros(2,2);
b=zeros(1,2);
A(2,1)=1;
b(1)=0.5;
b(2)=0.5;

u=rkex(A, b, u0, T, cfl);      
      


