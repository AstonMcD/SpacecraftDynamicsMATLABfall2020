%Earth and Mars in 2 Body Orbits (Keplerian) around the Sun.
%heliocentric Ecliptic J2000 frame. for all vectors and coordinates

%||INITIAL CONDITIONS||%
%Gravitiational Parameters
muSun = 1.327e11; %km^3/s^2
muMars = 4.282e4; %km3/s2
muEarth = 398600; %km3/s2
radMars=3396; %km
%Positon(km) and Velocity(km/s) of Celestial Bodies%
%time t=t0 on Jan 1st,2020
%   rv__ = [x, y, z, vx, vy, vz] (OR) [x,y,z,u,v,w]
rvEi = [-2.5527e+07, 1.4486e+08, -6.5836e+03, -2.9823e+01, -5.2819e+00, 3.1683e-04];
rvMi = [-1.9717e+08, -1.3290e+08, 2.0529e+06, 1.4449e+01, -1.8018e+01, -7.3209e-01];

%VARIABLE INDEX
%   i ~ initial
%   E ~ Earth
%   M ~ Mars
%   L ~ Launch
%   A ~ Arrival
%   sc ~ spacecraft
%   P ~ Plot/used in plot
%   I ~ Interpolated
%   MC ~ Mars Centric

%oe = [a e i w bigW v]
oe_Ei = hw6rv2oe(rvEi, muSun);
oe_Mi = hw6rv2oe(rvMi, muSun);

%Confirm Period of Planets is accurate
T_Esec = 2*pi*sqrt((oe_Ei(1))^3/muSun); 
T_Eday = T_Esec/(60*60*24); 
T_Msec = 2*pi*sqrt((oe_Mi(1))^3/muSun);
T_Mday = T_Msec/(60*60*24); 
%come out to be reasonably close, but not as accurate as I would hope.
%% Problem 1.
%Propogate Earth Orbit to time of launch (July 30, 2020 [t=t0+211 days])
timeLaunch = 211*86400; 
err=odeset('AbsTol',1e-8,'RelTol',1e-8);
[tE, rvE_1]=ode45(@(t,r) body2(t,r,muSun), [0 timeLaunch], rvEi, err);
[tM, rvM_1]=ode45(@(t,r) body2(t,r,muSun), [0 timeLaunch], rvMi, err);
rvEL=rvE_1(end,:);
rvML=rvM_1(end,:);
% EARTH: orbital elements with eccentric and true anomalies
oeEL=hw6rv2oe(rvEL,muSun);
EEL=2*atan(sqrt((1-oeEL(2))/(1+oeEL(2)))*tan(oeEL(6)/2));

%Propogate Mars and Earth Orbit to time of arrival (Feb 18, 2021 [t=t0+211+203 days])
timeArrival=414*86400;
[tE, rvE_2]=ode45(@(t,r) body2(t,r,muSun), [0 timeArrival], rvEi, err);
[tM, rvM_2]=ode45(@(t,r) body2(t,r,muSun), [0 timeArrival], rvMi, err);
rvEA=rvE_2(end,:);
rvMA=rvM_2(end,:);
% MARS: orbital elements with eccentric and true anomalies
oeMA=hw6rv2oe(rvMA,muSun);
EMA=2*atan(sqrt((1-oeMA(2))/(1+oeMA(2)))*tan(oeMA(6)/2));

%% PROBLEM 2.
%From Launch to Arrival is 203 days
timeFlight = 203*84600;
%initial guess w/ hohman 
%[I tried 2 different hohman functions one (with and without) capture]
[dv,tof,dva,dvp] = hohman2(muSun, norm(rvEL(1:3)), norm(rvMA(1:3))); %(mu, r1, r2)
%Vs/c and Vinf at launch
VinfL=dvp*rvEL(4:6)./norm(rvEL(4:6));
rvscL = [rvEL(1) rvEL(2) rvEL(3) VinfL(1) VinfL(2) VinfL(3)];
rvscLsave=rvscL;

%error in x,y, and z directions
errX=0; errY=0; errZ=0;
%Tolerance and Errors
tolUpdate=1e5;
errsc = 1e8;
tolFinal=1;
errVloop=1e7;
% Iteration, Trial&Error
while errsc>tolFinal && tolUpdate>1e-6
    %set up change in velocity(step)
    step = tolUpdate/1e7;         
    i=0;
while errVloop>tolUpdate && i<10000
    %change in velocity x-direction
    if errX > 1
        rvscL(4)=rvscL(4)+step;
    elseif errX < -1
        rvscL(4)=rvscL(4)-step;
    else
        rvscL(4)=rvscL(4);
    end
    %change in velocity y direction
    if errY > 1
        rvscL(5)=rvscL(5)+step;
    elseif errY < -1
        rvscL(5)=rvscL(5)-step;
    else
        rvscL(5)=rvscL(5);
    end
    %change in velocity z direction
    if errZ > 1
        rvscL(6)=rvscL(6)+step;
    elseif errZ < -1
        rvscL(6)=rvscL(6)-step;
    else
        rvscL(6)=rvscL(6);
    end
    %propogate s/c trajectory and cross check with calculated Mars Arrival
    [tsc,rvsc]=ode45(@(t,r) body2(t,r,muSun), [0 timeFlight], rvscL, err);
    rvscA=rvsc(end,:);
    errX=rvMA(1)-rvscA(1);
    errY=rvMA(2)-rvscA(2);
    errZ=rvMA(3)-rvscA(3);
    errVloop=norm(rvscA(1:3)-rvMA(1:3));
    i=i+1;
end
errsc=errVloop;
tolUpdate = tolUpdate/100;
end
%compare Mars position w/ s/c Arrival
HoM = [abs(rvMA(1)-rvscA(1)), abs(rvMA(2)-rvscA(2)), abs(rvMA(3)-rvscA(3))];
HitOrMiss = norm(HoM);

%% PROBLEM 3.
%Propogate Earth and Mars orbits from Launch to Arrival
[teP,rvEarthL2AP]=ode45(@(t,r) body2(t,r,muSun), [timeLaunch timeArrival], rvEL, err);
[tmP,rvMarsL2AP]=ode45(@(t,r) body2(t,r,muSun), [timeLaunch timeArrival], rvML, err);

%Top View of combined trajectories
figure(1)
plot3(rvsc(:,1),rvsc(:,2),rvsc(:,3),"g",rvEarthL2AP(:,1),rvEarthL2AP(:,2),rvEarthL2AP(:,3),"b",rvMarsL2AP(:,1),rvMarsL2AP(:,2),rvMarsL2AP(:,3),"r");
view(2)
title("TOP VIEW: Combined Resulting Trajectories of Bodies from Launch to Arrival")
xlabel("x-position (km)")
ylabel("y-position (km)")
zlabel("z-position (km)")
legend("Spacecraft","Earth","Mars")

%3D View of combined trajectories
figure(2)
plot3(rvsc(:,1),rvsc(:,2),rvsc(:,3),"g",rvEarthL2AP(:,1),rvEarthL2AP(:,2),rvEarthL2AP(:,3),"b",rvMarsL2AP(:,1),rvMarsL2AP(:,2),rvMarsL2AP(:,3),"r");
title("3D VIEW: Combined Resulting Trajectories of Bodies from Launch to Arrival")
xlabel("x-position (km)")
ylabel("y-position (km)")
zlabel("z-position (km)")
legend("Spacecraft","Earth","Mars")

%% PROBLEM 4.
%interpolate spacecraft to Mars Centric frame
tscIP = tsc+timeLaunch;
%---Used for Animation------------
rvMIP=interp1(tmP,rvMarsL2AP,teP);
rvscIP=interp1(tscIP,rvsc,teP);
%---------------------------------
rEarthMCP=rvEarthL2AP(:,1:3)-rvMIP(:,1:3);
rscMCP=rvscIP(:,1:3)-rvMIP(:,1:3);

%,rEarthMCP(:,1),rEarthMCP(:,2),rEarthMCP(:,3),"b"
%Top View of s/c in Mars Centric Frame
figure(3)
plot3(rscMCP(:,1),rscMCP(:,2),rscMCP(:,3),"g",rEarthMCP(:,1),rEarthMCP(:,2),rEarthMCP(:,3),"b",0,0,0,"ro");
view(2)
title("TOP VIEW:Spacecraft Trajectory in Mars Centric Frame with Earth for reference")
xlabel("x-position (km)")
ylabel("y-position (km)")
zlabel("z-position (km)")
legend("Spacecraft","Earth","Mars")

%3D View of s/c in Mars Centric Frame
figure(4)
plot3(rscMCP(:,1),rscMCP(:,2),rscMCP(:,3),"g",rEarthMCP(:,1),rEarthMCP(:,2),rEarthMCP(:,3),"b",0,0,0,"ro");
title("3D VIEW: Spacecraft Trajectory in Mars Centric Frame with Earth for reference")
xlabel("x-position (km)")
ylabel("y-position (km)")
zlabel("z-position (km)")
legend("Spacecraft","Earth","Mars")

%% PROBLEM 5.
%Create Vs/c, Vbody, and Vinf vectors at Launch
vEL=[0 0 0;rvEL(4:6)];
vscL=[0 0 0;rvscL(4:6)];
vinfLvec=vscL-vEL+[rvEL(4:6);rvEL(4:6)];
magVinfLaunch=norm(vscL(2,:)-vEL(2,:));

%Plot the Vinf Launch vector
figure(5)
plot3(vEL(:,1),vEL(:,2),vEL(:,3),"b",vscL(:,1),vscL(:,2),vscL(:,3),"g",vinfLvec(:,1),vinfLvec(:,2),vinfLvec(:,3),"r")
grid on
title("3D VIEW: Vinf Launch Velocity Diagram")
xlabel("x-position (km)")
ylabel("y-position (km)")
zlabel("z-position (km)")
legend("V-body (Earth)","V-s/c","V-inf")

%Create Vs/c, Vbody, and Vinf vectors at Arrival
vMA=[0 0 0;rvMA(4:6)];
vscA=[0 0 0;rvscA(4:6)];
vinfAvec=vscA-vMA+[rvMA(4:6);rvMA(4:6)];
magVinfArrival=norm(vscA(2,:)-vMA(2,:));

%Plot the Vinf Arrival vector
figure(6)
plot3(vMA(:,1),vMA(:,2),vMA(:,3),"b",vscA(:,1),vscA(:,2),vscA(:,3),"g",vinfAvec(:,1),vinfAvec(:,2),vinfAvec(:,3),"r")
grid on
title("3D VIEW: Vinf Arrival Velocity Diagram")
xlabel("x-position (km)")
ylabel("y-position (km)")
zlabel("z-position (km)")
legend("V-body (Mars)","V-s/c","V-inf")

%% PROBLEM 6.
%used code from mathworks forum
%Animation for 3D view
figure(7)
filename = '3dTrajectory.gif';
for j = 1:length(teP)
    plot3(rvscIP(1:j,1),rvscIP(1:j,2),rvscIP(1:j,3),"g",rvEarthL2AP(1:j,1),rvEarthL2AP(1:j,2),rvEarthL2AP(1:j,3),"b",rvMIP(1:j,1),rvMIP(1:j,2),rvMIP(1:j,3),"r");
    drawnow
    title("ANIMATION (3D): Combined Resulting Trajectories of Bodies from Launch to Arrival")
    xlabel("x-position (km)")
    ylabel("y-position (km)")
    zlabel("z-position (km)")
    legend("Spacecraft","Earth","Mars")
    frame = getframe(7);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if j==1
        imwrite(imind,cm,filename,'gif','Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
end

%Animation for Top View
figure(8)
filename = 'TopTrajectory.gif';
for j = 1:length(teP)
    plot3(rvscIP(1:j,1),rvscIP(1:j,2),rvscIP(1:j,3),"g",rvEarthL2AP(1:j,1),rvEarthL2AP(1:j,2),rvEarthL2AP(1:j,3),"b",rvMIP(1:j,1),rvMIP(1:j,2),rvMIP(1:j,3),"r");
    view(2)
    drawnow
    title("ANIMATION (Top): Combined Resulting Trajectories of Bodies from Launch to Arrival")
    xlabel("x-position (km)")
    ylabel("y-position (km)")
    zlabel("z-position (km)")
    legend("Spacecraft","Earth","Mars")
    
    frame = getframe(8);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if j==1
        imwrite(imind,cm,filename,'gif','Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
end

%% PROBLEM 7.
khat=[4.461321045940511e-01  -5.583636558975008e-02 8.932236257109472e-01];
    % h = r x v
h=cross(rvscA(1:3),rvscA(4:6));
    %inclination
isc=acosd(dot(h,khat)/(norm(h)*norm(khat)));

%% PROBLEM 8.
C3 = magVinfLaunch^2;
m0 = 4250;
vinfA=vscA(2,:)-vMA(2,:);
%collinear manuever at mars capture
vpM=sqrt(magVinfArrival^2 + 2*muMars/(radMars+1000));
vcM=sqrt(muMars/(radMars+1000));
deltaV=vpM-vcM;
g0=9.8065e-3;
lsp = 300; % seconds
vE=g0*lsp;
mL=m0*exp(-1*deltaV/vE);
mP=m0-mL;


