function y = expA2_vec(eign_expA2,v)
% Compute exp(-tau*A2)*v
% [n,m] = size(eign_expA);
% rex = reshape(v,Nx,Ny);
rex = fft2(v);
rex = rex.*eign_expA2;
y = real(ifft2(rex));
