function [U,S,V,r] = truncated_svd(A0,tol)
% [m,n] = size(A0);
% tol = [rtol;atol];
[U,sing_vals,V] = svd(A0,'vector');
rtol = tol(1); atol = tol(2);
if isempty(rtol)
    r = sum(sing_vals > atol);
else
    r = sum(sing_vals > sing_vals(1)*rtol);
end
if r == 0 % trivial case
    U = []; 
    S = [];
    V = [];
else % general case
    U = U(:,1:r);
    S = diag(sing_vals(1:r));
    V = V(:,1:r);
end