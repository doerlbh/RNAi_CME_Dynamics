% RNA interference modeled by Deterministic ODE
% Author: Baihan Lin
% Date: Nov 2016

%  dvdt = ep1 * (1 - v - H1 * v * r);
%  drdt = ep1 * (1 - r - H2 * v * r);
%  dqdt = ep2 * (- q + H3 * v * r);
%  dsdt = r - s;



function dy = iRNA_ODE(t,y,p)
 
   ep1 = p(1);
   ep2 = p(2);
   

   dy = zeros(4,1);
   
end