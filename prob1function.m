function xd = prob1(t, r)
mu = 1;

xd = zeros(6, 1);
r3 = (sqrt((r(1)^2 + r(2)^2 + r(3)^2)))^3;

xd(1) = r(4);
xd(2) = r(5);
xd(3) = r(6);
xd(4) = (-1*mu/r3)*r(1);
xd(5) = (-1*mu/r3)*r(2);
xd(6) = (-1*mu/r3)*r(3);

end