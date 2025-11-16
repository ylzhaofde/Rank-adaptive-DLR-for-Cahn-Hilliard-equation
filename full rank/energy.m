function y = energy(u,hx,hy,Ax,Ay)
E1=(-0.5*hx*hy*0.01)*trace((u*Ay+Ax*u)*u');
E2=(0.25*hx*hy)*norm(u.^2-1,"fro")^2;
y=E1+E2;