function err = fitellipse(start, lbx, lby,cx,cy)
a = start(1);
b = start(2);
theta = start(3);
% cx = start(4);
% cy = start(5);
eq = (((lbx-cx).*cos(theta)-(lby-cy).*sin(theta))./a).^2 + (((lbx-cx).*sin(theta)+(lby-cy).*cos(theta))./b).^2;
err = sum((1-eq).^2);