function u0 = u1fun(x,y,epsilon)
% initial value
Nx = length(x) - 1;
Ny = length(y) - 1;
u0 = zeros(Nx + 1,Ny + 1);
% star: (x,y) \in [0,1]^2, kappa = 1e-4, T = 0.1
% for j = 1:Ny + 1
%     for i = 1:Nx + 1
%         if x(i) > 0.5
%             theta = atan2(y(j) - 0.5, x(i) - 0.5);
%             u0(i,j) = tanh((0.25 + 0.1 * cos(6 * theta) - ...
%                         (sqrt((x(i) - 0.5)^2 + (y(j) - 0.5)^2))) / (sqrt(2) * epsilon));    
%         else
%             theta = pi + atan2(y(j) - 0.5, x(i) - 0.5);
%             u0(i,j) = tanh((0.25 + 0.1 * cos(6 * theta) - ...
%                         (sqrt((x(i) - 0.5)^2 + (y(j) - 0.5)^2))) / (sqrt(2) * epsilon));                   
%         end
%     end
% end

% torus: (x,y) \in [-1,1]^2, kappa = 1e-4, T = 0.1
% R1 = 0.4; R2 = 0.3;  
% for j = 1:Ny + 1
%     tmp = (x - 0).^2 + (y(j) - 0).^2;
%     u0(:,j) = -1 + tanh((R1 - sqrt(tmp))/(sqrt(2) * epsilon)) - ...
%                   tanh((R2 - sqrt(tmp))/(sqrt(2) * epsilon));
% end

% dumbbell: (x,y) \in [0,2] x [0,1], kappa = 1e-4, T = 0.1
% R0 = 0.2;
% for j = 1:Ny + 1
%     for i = 1:Nx + 1
%         if x(i) > 0.4 && x(i) < 1.6 && y(j) > 0.4 && y(j) < 0.6
%             u0(i,j) = 1;
%         else
%             u0(i,j) = 1 + tanh((R0 - sqrt((x(i) - 0.3).^2 + (y(j) - 0.5).^2))/(sqrt(2) * epsilon)) + ...
%                         tanh((R0 - sqrt((x(i) - 1.7).^2 + (y(j) - 0.5).^2))/(sqrt(2) * epsilon));
%         end
%     end
% end
%maze: (x,y) \in [0,1]^2, kappa = 1e-4, T = 0.1
for i = 1:Nx + 1
    for j = 1: Ny + 1
        u0(i, j) = -1;
        if x(i) > 0.1 && x(i) < 0.2 && y(j) > 0.1 && y(j) < 0.9 || ...
           x(i) > 0.2 && x(i) < 0.8 && y(j) > 0.1 && y(j) < 0.2 || ...
           x(i) > 0.8 && x(i) < 0.9 && y(j) > 0.1 && y(j) < 0.8 || ...
           x(i) > 0.3 && x(i) < 0.8 && y(j) > 0.7 && y(j) < 0.8 || ...
           x(i) > 0.3 && x(i) < 0.4 && y(j) > 0.3 && y(j) < 0.7 || ...
           x(i) > 0.4 && x(i) < 0.7 && y(j) > 0.3 && y(j) < 0.4 || ...
           x(i) > 0.6 && x(i) < 0.7 && y(j) > 0.4 && y(j) < 0.6
        elseif x(i) > 0.1 && x(i) < 0.15 && y(j) > 0.1 && y(j) < 0.95 || ...
           x(i) > 0.1 && x(i) < 0.9 && y(j) > 0.1 && y(j) < 0.15 || ...
           x(i) > 0.85 && x(i) < 0.9 && y(j) > 0.1 && y(j) < 0.8 || ...
           x(i) > 0.3 && x(i) < 0.9 && y(j) > 0.75 && y(j) < 0.8 || ...
           x(i) > 0.3 && x(i) < 0.35 && y(j) > 0.3 && y(j) < 0.8 || ...
           x(i) > 0.3 && x(i) < 0.7 && y(j) > 0.3 && y(j) < 0.35 || ...
           x(i) > 0.65 && x(i) < 0.7 && y(j) > 0.3 && y(j) < 0.6 || ...
           x(i) > 0.45 && x(i) < 0.7 && y(j) > 0.55 && y(j) < 0.6 || ...
           x(i) > 0.45 && x(i) < 0.5 && y(j) > 0.4 && y(j) < 0.6 || ...
           x(i) > 0.45 && x(i) < 0.6 && y(j) > 0.4 && y(j) < 0.45
             u0(i, j) = 1;
        end
    end
end
