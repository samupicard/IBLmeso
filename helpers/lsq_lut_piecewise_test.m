% Test for lsq_lut_piecewise function

% generate noisy x and y data
x = [-2:0.1:2]';
x = x+2*(rand(size(x))-0.5);
x = max(min(x,3),-3);

y = x+x.^2+2*sin(x)+5*(rand(size(x))-0.5);

% vector of 1-D look-up table "x" points
XI = linspace(min(x),max(x),3);

% obtain vector of 1-D look-up table "y" points
YI = lsq_lut_piecewise( x, y, XI );

% plot fit
plot(x,y,'.',XI,YI,'+-')
legend('experimental data (x,y(x))','LUT points (XI,YI)')
title('Piecewise 1-D look-up table least square estimation')
