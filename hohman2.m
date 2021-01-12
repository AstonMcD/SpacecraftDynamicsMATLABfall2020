function [dv,tof,va,vp] = hohman2(mu,r1,r2)
vt1 = sqrt(((-2*mu)/(r1+r2))+(2*mu/r1));
vc1=sqrt(mu/r1);
vt2 = sqrt(((-2*mu)/(r1+r2))+(2*mu/r2));
vc2=sqrt(mu/r2);

dvp = abs(vt1-vc1);
vp=sqrt(dvp^2 + (2*398600/6378));
dva = abs(vt2-vc2);
va = sqrt(dva^2 + (2*42828/3396));
dv = va + vp;
tof = pi*sqrt(((r1+r2)^3)/(8*mu));
end