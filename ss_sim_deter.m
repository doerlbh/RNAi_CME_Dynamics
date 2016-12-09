% RNA interference modeled by Deterministic ODE
% Author: Baihan Lin
% Date: Nov 2016

clear all; close all

gm1 = 
gm2 = 
gm3 = 
K1 = 
K2 =
H = 

ep1 = gm2/gm3;
ep2 = gm1/gm3;
H1 = H*K1/(gm1*gm3);
H2 = H*K2/(gm1*gm3);
H3 = H^2*ep1*K1*K2/(gm2^2);

for numreps=1:1

p = [ep1, ep2, H1, H2, H3];
x0 = [0, 0, 0, 0];

Tmax=80;
options = odeset('NonNegative',1:4); %make solutions nonnegative

[T,Y] = ode45(@(t,y) iRNA_ODE(t,y,p),[0 Tmax],x0,options);

figure(1)
set(gca,'FontSize',16)
3Dplot(T,Y(:,1:2),'LineWidth',3); hold on;
3Dplot(T,Y(:,4),'LineWidth',3); hold on;
legend('v (iRNA)','r (mRNA)','s (protein)')
xlabel('t'); 
ylabel('%');

end