function x = StateA(x,u)

Fyf = x(3);
Fyr = x(4);
Fxf = x(5);
beta = u(1);
deltaf = u(2);
m = 1600;
izz = 2310;
lf = 1.11;
lr = 1.756;
vx_p =  (Fxf*cos(beta - deltaf) + Fyf*sin(beta - deltaf) + Fyr*sin(beta))/m;
yawrate_p = ((Fyf*cos(deltaf) + Fxf*sin(deltaf))*lf - Fyr*lr)/izz;
x(1) = x(1) + vx_p * 0.01;
x(2) = x(2) + yawrate_p * 0.01;
x(3) = x(3);
x(4) = x(4);
x(5) = x(5);
end