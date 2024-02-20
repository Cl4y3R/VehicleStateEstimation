function y = MeasureA(x,u)
m = 1600;
y = zeros(4,1);
y(1) = x(1);
y(2) = x(2);
Fyf = x(3);
Fyr = x(4);
Fxf = x(5);
deltaf = u(2);
y(3) = (-Fyf*sin(deltaf) + Fxf*cos(deltaf))/m;
y(4) = (Fyf*cos(deltaf) + Fyr + Fxf*sin(deltaf))/m;
end