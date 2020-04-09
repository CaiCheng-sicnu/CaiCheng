%Graphene tight binding (TB) model for 3D bands in matlab.
% % % %plot 3D bands
x=linspace(-2.9,2.9,101); 
y=linspace(-12,12,151); 
[a,b]=meshgrid(x,y);
z= 3.033*sqrt(1+4*cos(b*0.1*pi).*cos(a*1.2283)+4*(cos(a*1.2283)).^2); 
zz=-3.033*sqrt(1+4*cos(b*0.1*pi).*cos(a*1.2283)+4*(cos(a*1.2283)).^2); 
figure 
surf(a,b,z) 
title('Energy dispersion relation of graphene') 
xlabel('K(1/A)'),ylabel(' Energy of graphene(ev)') 
hold on 
surf(a,b,zz)
