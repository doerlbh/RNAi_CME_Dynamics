% RNA interference modeled by Monte Carlos Simulation
% Author: Baihan Lin
% Date: Nov 2016

clear all; close all

p1=0.0737;  %k25
p2=1.29;    %k32
p3=0.0411;  %k72
p4=0.0212;  %k43
p5=0.23;    %k64
p6=0.101;   %k56
p7=0.23;    %k57

Nreceptor = 100000;
p = [p1; p2; p3; p4; p5; p6; p7];
x0 = [0; Nreceptor; 0; 0; 0; 0; 0];

rng(1);

Nstep=10000000;
stptime = zeros(Nstep,1);
time = zeros(Nstep,1);
xall = zeros(Nstep,7);
xall(1,:) = x0;
x = x0;

for step = 1 : Nstep - 1
   [xnew, tau] = new_stochastic_update(x, p);
   x = xnew;
   stptime(step+1) = tau;
   time(step+1) = time(step) + tau;
   xall(step+1,:) = x;
end

figure(1)
set(gca,'FontSize',16);
plot(time,100*xall(:,1:7)/Nreceptor,'LineWidth',3);
legend('x1','x2','x3','x4','x5','x6','x7');
title('Monte Carlos Stochastic Dynamics of 10000 Receptors vs. time')
xlabel('time (min)'); 
ylabel('% of receptors');
