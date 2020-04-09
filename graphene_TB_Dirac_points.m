%Graphene tight binding (TB) model for K and K' Dirac point.
%%% Hamiltion
syms qx qy M
M=0.001;
%h=[M qx-1i*qy; qx+1i*qy -M]
%eig(h);
%   (M.^2 + qx.^2 + q.y^2).^(1/2)
%  -(M.^2 + qx.^2 + qy.^2).^(1/2)
% % % %plot 3D bands
x=linspace(-0.2,0.2,101); 
y=linspace(-0.2,0.2,101); 
[qx,qy]=meshgrid(x,y);
z= (M.^2 + qx.^2 + qy.^2).^(1/2); 
zz=-(M.^2 + qx.^2 + qy.^2).^(1/2); 
figure 
surf(qx,qy,z) 
hold on 
surf(qx,qy,zz)
hold on
title('Dirac point of graphene') 
xlabel('Kx(1/A)'),ylabel('Ky(1/A)') 
zlabel('Energy (ev)') 
