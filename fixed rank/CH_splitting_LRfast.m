% fixed low-rank splitting method for CH equation
% Example 1 in the article(in Chinese):
% S. Wang, Y.-L. Zhao, X.-M. Gu, 
%  A RANK-ADAPTIVE LOW-RANK SPLITTING METHOD FOR THE CAHN-HILLIARD EQUATION

clear; close all;
clc;

x_L = 0; x_R = 2*pi;
y_L = 0; y_R = 2*pi;
T = 0.1; % time
kappa = 0.01; % for Example 1; for Example 2, kappa = 0.1
M = 2^9; % time steps
Nx = 2^8; % X-direction
Ny = Nx;  % Y-direction
fu = @(u) u.^3 - u;
approx_r = 2;  % rank

hx = (x_R - x_L)/Nx;
hy = (y_R - y_L)/Ny;
tau = T/M;
x = x_L + (0:Nx)'*hx;
y = y_L + (0:Ny)'*hy;
t = (0:M)'*tau;

u = zeros(Nx + 1,Ny + 1, M + 1);
% initial value
%example 2
%u(:,:,1)=u1fun(x,y,kappa);
%example 1
for j = 1:Ny + 1
       u(:,j,1) = 0.05*sin(x)*sin(y(j));
end
tic;
% the matrices
e1 = ones(Nx,1);
Ax = 1/hx^2*spdiags([e1, -2*e1, e1], [-1,0,1], Nx,Nx);
Ax(1,end) = 1/hx^2; Ax(end,1) = 1/hx^2;
Ay = 1/hy^2*spdiags([e1, -2*e1, e1], [-1,0,1], Ny,Ny);
Ay(1,end) = 1/hy^2; Ay(end,1) = 1/hy^2;
% The matrices
eign_Lapx = 1/hx^2*[0;-4*sin((1:Nx - 1)'*pi/Nx).^2];
eign_Lapy = 1/hy^2*[0;-4*sin((1:Ny - 1)'*pi/Ny).^2];
% eign_B = (1/hx^2*eign_Bx*ones(1,Ny) + 1/hy^2*ones(Nx,1)*eign_By.');
eign_Lap = eign_Lapx*ones(1,Ny) + ones(Nx,1)*eign_Lapy.';
eign_expA2 = exp(-tau*kappa*(eign_Lap).^2);
expA2vec = @(v) expA2_vec(eign_expA2,v);

% lowr-rank approx for initial value
[U,S,V] = svds(u(1:Nx,1:Ny,1),approx_r);
u(1:Nx,1:Ny,1) = U*S*V';
for j = 1:M
    % linear part:
    v = expA2vec(U*S*V');
    [U_v,S_v,V_v] = DLRA_proj_splitting(U,S,V,v);
    % nonlinear part:
    v = U_v*S_v*V_v';
    [U,S,V] = DLRA_proj_splitting(U_v,S_v,V_v,tau*(fu(v)*Ay + Ax*fu(v)) + v);
    u(1:Nx,1:Ny,j + 1) = U*S*V';
    % energy:
    en(j)=energy(u(1:Nx,1:Ny,j + 1),hx,hy,Ax,Ay);
end
u(end,:,2:end) = u(1,:,2:end);
u(:,end,2:end) = u(:,1,2:end);
interval = toc;

fprintf('Copyright (c) 2025 Y.-L. Zhao (School of Mathematical Sciences, Sichuan Normal University).\n');
fprintf('All rights reserved.\n');
fprintf('Fixed low-rank splitting for CH, Example 1\n');
fprintf('CPU: %.3f\n',interval);
