function [dv,tof,dva,dvp] = hohman(mu,r1,r2)
vt1 = sqrt(((-2*mu)/(r1+r2))+(2*mu/r1));
vc1=sqrt(mu/r1);
vt2 = sqrt(((-2*mu)/(r1+r2))+(2*mu/r2));
vc2=sqrt(mu/r2);

dvp = abs(vt1-vc1);
dva = abs(vt2-vc2);
dv = dva + dvp;
tof = pi*sqrt(((r1+r2)^3)/(8*mu));
end