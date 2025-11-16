function [U_1,S_1,V_1] = DLRA_proj_splitting(U_0,S_0,V_0,deltaA)
    
% Lie-Trotter splitting
% First order projector splitting integrator (KSL)
% for dynamical low-rank approximation 
% of the solution of a matrix differential equation 
% see PhD thesis: On dynamical low-rank integrators for matrix differential
% equations, Pages 33~34 for KSL with exact solution
%
% INPUT:
%      U,S,V: initial values for the factors 
%      tspan: time interval
% OUTPUT:
%      U,S,V: factors of the solution after one step

%%  KSL with exact solution

% temp = U_0*S_0*V_0';
% deltaA = temp./sqrt(exp(-2*tau) + (1 - exp(-2*tau))*temp.^2);
% K step: Update U
US1 = deltaA*V_0; % U*S
[U_1,~] = qr(US1,'econ');
% L step: Update V
VSH1 = deltaA'*U_1; 
[V_1,S_1] = qr(VSH1,'econ');
S_1 = S_1';

