function y = MeasureB(x,u)
y = zeros(3,1);
vx = u(1);
yawrate = u(2);
deltaf = u(3);
Fxf = u(4);
beta = x(1);
m = 1600;
lf = 1.11;
lr = 1.756;
cf = 100000;
cr = 100000;
y(1) = (cf + x(2))*(deltaf - beta - lf*yawrate/vx);
y(2) = (cr + x(3))*(- beta + lr*yawrate/vx);
y(3) = ((cf + x(2))*(deltaf - beta - lf*yawrate/vx)*cos(deltaf) + (cr + x(3))*(- beta + lr*yawrate/vx) + Fxf*sin(deltaf))/m;
end