% Graphene tight binding (TB) model for bands in matlab.
% matlab code in the following path: G-K-M-G is:
clear;clc
%% % you should input it
%cell vector in real space
syms a t 
a=2.44; % Angstrom
t=-2.5;
a1=[3*a/2, sqrt(3)*a/2];
a2=[3*a/2,-sqrt(3)*a/2];
delta1=[a/2, sqrt(3)*a/2];
delta2=[a/2,-sqrt(3)*a/2];
delta3=[-a,0];
K_point= [2*pi/3/a, 2*pi/3/sqrt(3)/a];
Kp_point=[2*pi/3/a,-2*pi/3/sqrt(3)/a];
b1=[2*pi/3/a, 2*sqrt(3)*pi/3/a];
b2=[2*pi/3/a,-2*sqrt(3)*pi/3/a];
Nt=101;

% plot band G-M (X1)
kx=linspace(0,2*pi/3/a,Nt);ky=0;
GM=[0,0;0,2*pi/3/a];
dist_GM=pdist(GM,'euclidean');
X1=linspace(0,dist_GM,Nt);
for Nk=1:Nt
    k=[kx(1,Nk) ky];
h12=(exp(1i*dot(k,delta1))+exp(1i*dot(k,delta2))+exp(1i*dot(k,delta3)));
h21=conj(h12);
H=-t*[0 h12;h21 0];
    [V,D]=eig(H);
    eigst=sum(D);
    E1(Nk,:)=sort(real(eigst));
end
plot(X1,E1,'b');
hold on;

% plot band M-K (X2)
kx=2*pi/3/a;ky=linspace(0,2*pi/3/sqrt(3)/a,Nt);
MK=[0,0;0,2*pi/3/a];
dist_MK=pdist(MK,'euclidean');
X2=linspace(0,dist_MK,Nt);
for Nk=1:Nt
    k=[kx ky(1,Nk)];
h12=(exp(1i*dot(k,delta1))+exp(1i*dot(k,delta2))+exp(1i*dot(k,delta3)));
h21=conj(h12);
H=-t*[0 h12;h21 0];
    [V,D]=eig(H);
    eigst=sum(D);
    E2(Nk,:)=sort(real(eigst));
end
plot(dist_GM+X2,E2,'b');
hold on;

% plot band K-G (X3)
kx=linspace(2*pi/3/a,0,Nt);ky=linspace(2*pi/3/sqrt(3)/a,0,Nt);
KG=[0,0;0,2*pi/3/a];
dist_KG=pdist(KG,'euclidean');
X3=linspace(0,dist_KG,Nt);
for Nk=1:Nt
    k=[kx(1,Nk) ky(1,Nk)];
h12=(exp(1i*dot(k,delta1))+exp(1i*dot(k,delta2))+exp(1i*dot(k,delta3)));
h21=conj(h12);
H=-t*[0 h12;h21 0];
    [V,D]=eig(H);
    eigst=sum(D);
    E3(Nk,:)=sort(real(eigst));
end
plot(dist_GM+dist_MK+X3,E3,'b');
hold on;

% plot the vertical line,
line([dist_GM dist_GM],[8,-8],'Color','k');
line([dist_GM+dist_MK dist_GM+dist_MK],[8,-8],'Color','k');
line([dist_GM+dist_MK+dist_KG dist_GM+dist_MK+dist_KG],[8,-8],'Color','k');
line([0,dist_GM+dist_MK+dist_KG],[0,0],'linestyle','--');
% set xlabel ticks
set(gca,'XTick',[0,dist_GM,dist_GM+dist_MK,dist_GM+dist_MK+dist_KG])
set(gca,'XTickLabel',{'\Gamma','M','K','\Gamma'});
%grid minor;
axis([0, dist_GM+dist_MK+dist_KG, -8, 8]);
xlabel('K(1/A)','Fontsize',14);
ylabel('Energy (eV)','Fontsize',14);
title('graphene TB','Fontsize',14);

