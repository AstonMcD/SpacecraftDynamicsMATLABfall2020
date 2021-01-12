function xd = prob0(t, x)
k=2;
m=3;

xd= zeros(2,1);
xd(1) = x(2);
xd(2) = -(k/m)*x(1);

end