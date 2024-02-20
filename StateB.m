function x = StateB(x,u)
beta = x(1);
vx = u(1);
yawrate = u(2);
deltaf = u(3);
Fxf = u(4);
m = 1600;
cf = 100000;
cr = 100000;
lf = 1.11;
lr = 1.756;
Fyf_est = (cf + x(2))*(deltaf - beta - lf*yawrate/vx);
Fyr_est = (cr + x(3))*(- beta + lr*yawrate/vx);
beta_p =  -yawrate + (Fxf*cos(deltaf - beta) + Fyf_est*cos(deltaf - beta) + Fyr_est*cos(beta))/(m*vx);
x(1) = x(1) + beta_p * 0.01;
x(2) = x(2);
x(3) = x(3);
end