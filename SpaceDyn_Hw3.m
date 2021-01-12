%PROBLEM 0%
k=2;
m=3;

period = abs(2*pi*(sqrt(m/k)));
x0 = [1 0];
r_init = 1;
t0 = 0;
tf = 5*period;
tspan = [t0 tf];
fprintf("prob0-i: Period rounded to 15 decimal places is %.14f \n", period)

[t, xd] = ode45(@(t,xd) prob0(t,xd), tspan, x0);
xd = xd(:,1);

figure(1);
plot(t, xd,"-o")
title("prob0-ii. Time vs Oscillator Position")
xlabel("Time (TU)")
ylabel("Position (LU)")

error0iii=xd(end,1)-r_init;
fprintf('prob0-iii. The error in r after 5 periods with a tolerance of 1e-3 is %.5f\n',error0iii)

err0 = zeros(1,13);
tol0 = zeros(1,13);

for i = 2:14
    tol = 1*10^(-i);
    options = odeset('RelTol', tol, 'AbsTol', tol);
    [t,xd] = ode45(@prob0, tspan, x0, options);
    err = abs((1-xd(length(xd))));
    
    tol0(i-1) = tol;
    err0(i-1) = err;
    
end

figure(2);
loglog(tol0, err0, 'r-o')
title("prob0-iv: Oscillator Error Tolerance vs Error Value")
xlabel("Error Tolerance")
ylabel("Error Value (LU)")

%PROBLEM 1%
%1 ORBIT%
mu = 1;
init1 = [1 1 1 0.7 0 0];
T = 5.003;
tspan1 = [0 T];
options1 = odeset('RelTol', 1e-14, 'AbsTol', 1e-14);

figure(3);
[t, xd] = ode45(@(t, r) prob1(t, r), tspan1, init1, options1);
plot3(xd(:,1), xd(:,2), xd(:,3), 'b-o')
axis equal
view(-40, 30)
title("prob1-i: Two Body Tragectory")
xlabel("x-Pos(LU)")
ylabel("y-Pos(LU)")
zlabel("z-Pos(LU)")
xd(end,1)
period1dt=linspace(0,T,100);
[t, xd] = ode45(@(t, r) prob1(t, r), period1dt, init1, options1);

figure(4);
plot3(xd(:,1), xd(:,2), xd(:,3), 'b-o')
axis equal
view(-40, 30)
title("prob1-ii: Two Body Tragectory + Time Spacing")
xlabel("x-Pos(LU)")
ylabel("y-Pos(LU)")
zlabel("z-Pos(LU)")

%ERROR%
error = abs((1-xd(length(xd))));
fprintf("prob1-iii: The Error after 1 period is %.14f \n", error)

%15 ORBITS%
period1dt15=linspace(0,15*T,1500);
[t, xd] = ode45(@(t, r) prob1(t, r), period1dt15, init1, options1);

figure(5);
plot3(xd(:,1), xd(:,2), xd(:,3), 'r-o')
axis equal
view(-40, 30)
title("prob1-iv: Two Body Tragectory for 15 periods + Time Spacing")
xlabel("x-Pos(LU)")
ylabel("y-Pos(LU)")
zlabel("z-Pos(LU)")

error1iv=sqrt((xd(end,1)-1)^2+(xd(end,3)-1)^2+(xd(end,5)-1)^2);
fprintf('prob1-iv: The error after fifteen periods is %.4f\n',error1iv)

%Error%
errInit1 = zeros(1,13);
tolInit1 = zeros(1,13);

for i = 2:14
    tol1 = 1*(10^(-i));
    options = odeset('RelTol', tol1, 'AbsTol', tol1);
    [t,r] = ode45(@(t, r) prob1(t, r), tspan1p4, init1, options);
    err1 = abs((1-r(length(r))));
    
    tolInit1(i-1) = tol1;
    errInit1(i-1) = err1;
    
end


figure(6);
loglog(tolInit1, errInit1, 'g-o')
xlabel("Error Tolerance")
ylabel("Error Value (LU)")
title("prob1-v: 2 Body System Error Tolerance vs Error Value ")

%PROBLEM 2%
G=1;
M1 = 1;
M2 = 2;
M3 = 1;
init2 = [1 0 0 0 0 0, 2 1 0.2 0 0 0, 3 0 0 0 2 0];
tspan2 = [0 5];

options2 = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
[t, r] = ode45(@(t,r) prob2function(t,r,G,M1,M2,M3), tspan2, init2, options2);

%CENTER OF MASS%
l=length(r);
sumMass=M1+M2+M3;
rcm=zeros(l,3);
for i= 1:l
    rcm(i,1:3)=(M1*r(i,1:3)+M2*r(i,7:9)+M3*r(i,13:15))/sumMass;
end
figure(7);
plot(t,rcm(:,1),'r',t,rcm(:,2),'b',t,rcm(:,3),'g')
title("prob2-i: Center of Mass Components")
legend("x-dir","y-dir","z-dir")
xlabel("Time (TU)")
ylabel("Position of CoM (LU)")

%TOTAL ANGULAR MOMENTIM%
H1=zeros(l,3);
for i= 1:l
    H1(i,1:3)=cross(r(i,1:3),M1*r(i,4:6));
end
H2=zeros(l,3);
for i= 1:l
    H2(i,1:3)=cross(r(i,7:9),M2*r(i,10:12));
end
H3=zeros(l,3);
for i= 1:l
    H3(i,1:3)=cross(r(i,13:15),M3*r(i,16:18));
end
hTot=H1+H2+H3;

figure(8);
plot(t,hTot(:,1),'r',t,hTot(:,2),'b',t,hTot(:,3),'g')
title("prob2-ii: Angular Momentum Components")
legend("x-dir","y-dir","z-dir")
xlabel("Time (TU)")
ylabel("Angular Momentum (MU*LU^2/TU)")

hTotxMean=hTot(:,1)-mean(hTot(:,1));
hTotyMean=hTot(:,2)-mean(hTot(:,2));
hTotzMean=hTot(:,3)-mean(hTot(:,3));
figure(9);
plot(t,hTotxMean,'b',t,hTotyMean,'r',t,hTotzMean,'g')
title("prob2-iii: Angular Momentum minus Average Components")
legend("x-dir","y-dir","z-dir")
xlabel("Time (TU)")
ylabel("Angular Momentum (MU*LU^2/TU)")

%ENERGY%
r12E=sqrt((r(:,1)-r(:,7)).^2 + (r(:,2)-r(:,8)).^2 + (r(:,3)-r(:,9)).^2);
r13E=sqrt((r(:,1)-r(:,13)).^2 + (r(:,2)-r(:,14)).^2 + (r(:,3)-r(:,15)).^2);
r23E=sqrt((r(:,13)-r(:,7)).^2 + (r(:,14)-r(:,8)).^2 + (r(:,15)-r(:,9)).^2);
K=0.5*M1*(r(:,4).^2+r(:,5).^2+r(:,6).^2) + 0.5*M2*(r(:,10).^2+r(:,11).^2+r(:,12).^2) + 0.5*M3*(r(:,16).^2+r(:,17).^2+r(:,18).^2);
U=G*M1*M2./r12E+G*M1*M3./r13E+G*M3*M2./r23E;

E=K-U;
figure(10);
plot(t,K,'r',t,U,'b',t,E-mean(E),'k')
title("prob2-iv: Energy Components")
legend("KE","PE","Total E - mean(Total E)")
xlabel("Time (TU)")
ylabel("Energy (MU*LU^2/TU^2)")

%Trajectories%
figure(11);
plot3(r(:,1),r(:,2),r(:,3),'k',r(:,7),r(:,8),r(:,9),'r',r(:,13),r(:,14),r(:,15),'g',rcm(:,1),rcm(:,2),rcm(:,3),'c')
axis equal
view(0,90)
title("prob2-v: Three Body Trajectories & Center of Mass")
legend("Body-1","Body-2","Body-3","CoM")
xlabel("x-Pos(LU)")
ylabel("y-Pos(LU)")
zlabel("z-Pos(LU)")

r21_T=r(:,7:9)-r(:,1:3);
figure(12);
plot3(r21_T(:,1),r21_T(:,2),r21_T(:,3))
axis equal
view(-40,30)
title("prob2-vi: Trajectory of body 2 with respect to body 1")
xlabel("Relative x-Pos(LU)")
ylabel("Relative y-Pos(LU)")
zlabel("Relative z-Pos(LU)")
