function [U_1,S_1,V_1,rmax] = adapDLRA_proj_splitting(U_0,S_0,V_0,deltaA)
    
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

% % K step: Update U
% US1 = deltaA*V_0; % U*S
% [U_1,~] = qr([US1,U_0],0);
% % L step: Update V
% VSH1 = deltaA'*U_0; % V*S^H 
% % or use the following command for smaller relerr, but cannot parallel
% % VSH1 = deltaA'*U_1; 
% [V_1,~] = qr([VSH1,V_0],0);
spmd
    if spmdIndex == 1 % Worker 1: K step (Update U)
        US1 = deltaA * V_0; % U*S
        [U_1, ~] = qr([US1, U_0], 'econ');
    elseif spmdIndex == 2 % Worker 2: L step (Update V)
        VSH1 = deltaA' * U_0; % V*S^H
        [V_1, ~] = qr([VSH1, V_0], 'econ');
    end
end
% Collect results from workers
U_1 = U_1{1}; % Extract U_1 from worker 1
V_1 = V_1{2}; % Extract V_1 from worker 2
% S step: Update S
S_1 = U_1'*deltaA*V_1;

    
    % Compute singular values of S1 and decide how to truncate:
    [U,S,V] = svd(S_1);
    
    tol = 1e-2;
%     tol = 1e-4*norm(S);
    
    rmax = size(S_1,1)/2;
    sg = diag(S); %sum(sg)
    
    for j=1:2*rmax
        tmp = sqrt(sum(sg(j:2*rmax)).^2);
        if(tmp<tol)
            break;
        end
    end
    
    rmax = j;
    rmax = min(rmax,30);
    
    % To use in the rank-fixed way - just comment the previous part.
    
    %Truncation:
    U_1 = U_1*U;
    V_1 = V_1*V;
    
    S_1 = S(1:rmax,1:rmax);
    U_1 = U_1(:,1:rmax);
    V_1 = V_1(:,1:rmax);  

end


