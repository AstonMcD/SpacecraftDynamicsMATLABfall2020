%% Solving the orbit
clear all

tMax = 603148; % seconds in 7 revolutions

% M = body revolutions
% N = orbits
% Seconds in one revolution is 86164
% Maximum M value = (604800 * (2pi / 86164)) / 2pi = 7.0191724 

muEarth = 3.986e5;
wB = (2*pi)/86164;
M = 7;
N = 1:200;

%placeholder to store correct
a0 = [0];

% semi major axis must 6678 <= a =< 6683.
for i = 1:length(N)
   
sma(i) = muEarth^(1/3) * (M/(N(i)*wB))^(2/3);
    
    if sma(i) >= 6678
        if sma(i) <= 6683
            Norbits = i;   
        a0(end+1) = sma(i);
        end
    end
          
end


  a = a0(2);

% I Set M to a value under 7.0191724 which is the exact revolutions in 1 week. 
% range for semi major axis (a). 6678 <= a => 6683.

Trgt = (2*pi) * sqrt((a^3)/muEarth) % Time of period in seconds

Trepeat = Norbits * Trgt % Time of repeat period in seconds

%% INSERT GROUNDTRACK TOOL CODE

toolInput = [a 0 85 0 0 0, Norbits 0 0 0];
%           [a e i omega LAN nu0, N thetaG0 ObsLon ObsLat]

%----------------------------Tool Input Index-----------------------------%
    a = toolInput(1);               %semi-major axis (km)
    e = toolInput(2);               %eccentricity
    i = deg2rad(toolInput(3));      %inclination (deg)
    omega = deg2rad(toolInput(4));  %argument of periapse (deg)
    LAN = deg2rad(toolInput(5));    %longitude of ascending node (deg)
    nu0 = deg2rad(toolInput(6));    %initial true anomaly(deg)
    N = toolInput(7);               %number of complete spacecraft orbits to propagate
    thetaG0 = toolInput(8);         %greenwich sidereal angle at t=t0 (deg)
    Lon = toolInput(9);             %observer longitude (deg)
    Lat = toolInput(10);            %observer latitude (deg)
%-------------------------------------------------------------------------%
muEarth = 3.986e5; %(km^3/s^2)
radEarth = 6378; %(km)
wB = (2*pi)/86164; %Earth Rotation

%% PROBLEM 1. Ground Track on Longitude Latitude Map
%Period of orbit
T = 2*pi * sqrt(a^3/muEarth); 
%Time to orbit
Ttotal= N * T;
tspan = 0:60:Ttotal; %minutes
%Set Position and Velocity
oe = [a e i omega LAN nu0];
rvSet = hw6oe2rv(oe,muEarth);

%Propogate orbit
err=odeset('AbsTol',1e-8,'RelTol',1e-8);
[t, rv]=ode45(@(t,r) body2(t,r,muEarth), tspan, rvSet, err);
rvFinal=rv(end,:);
%Calculate GST position
thetaGST = thetaG0 + (wB * tspan);

% tol = 1e-12;
% options = odeset('Reltol',tol,'Abstol',tol);
% 
% [t,rv] = ode45(@TwoBodyODE, tspan, rv0, options);

% Rotating the Frame and calculating Latitude and Longitude
for n = 1:length(rv)
    
    R3 = [cos(thetaGST(n)) sin(thetaGST(n)) 0; 
        -sin(thetaGST(n)) cos(thetaGST(n)) 0; 
        0 0 1];
    
    rECI = [rv(n,1);rv(n,2);rv(n,3)];
    
    rBody(n,1:3) = (R3*rECI);
    rNorm = norm(rBody(n,1:3));
    vecLat(n) = asin((rBody(n,3))/rNorm);
    vecLon(n) = atan2(rBody(n,2),rBody(n,1));
    
end

%convert back to degrees
vecLat = rad2deg(vecLat);
vecLon = rad2deg(vecLon);

%Set Image
img = imread('earthmap.jpg'); 
display = image(xlim,-[-90 90],img); 
uistack(display,'bottom')
axis([-180 180 -90 90])


%% PROBLEM 2. Body-Fixed ECEF Trajectory w/ Scaled Earth

%NOW WITH TOP VIEW
figure(1)
%Set Earth Image to Sphere
[x2,y2,z2] = sphere(40);
Earth = surf(x2 * radEarth, y2 * radEarth, z2 * radEarth);
img2 = flip(img,1);
set(Earth,'CData', img2,'FaceColor','texture');
hold on

%Plot 3D trajectory
plot3(rv(:,1), rv(:,2), rv(:,3),'-r');
view(2)
axis equal
title('Body-Fixed Frame Trajectory')
xlabel('X km')
ylabel('Y km')
zlabel('Z km')
legend('Earth','Trajectory')
hold off

%pause at the end of each section to give plots time to load
pause(1)
%% PROBLEM 3. Inertial Trajectory w/ Scaled Earth

figure(2)
%Set Earth Image to Sphere
[x3,y3,z3] = sphere(40);
Earth = surf(x3 * radEarth, y3 * radEarth, z3 * radEarth);
img3 = flip(img,1);
set(Earth,'CData', img3,'FaceColor','texture');
hold on

%Plot 3D trajectory
plot3(rBody(:,1), rBody(:,2), rBody(:,3),'-r');
axis equal
title('Inertial Frame Spacecraft Trajectory')
xlabel('X km')
ylabel('Y km')
zlabel('Z km')
legend('Earth','Trajectory')
hold off

%pause at the end of each section to give plots time to load
pause(1)

%% Why this trajectory?

% A trajectory of this design would be most useful for taking a scan of the
% entire globe perhaps to produce a map or a program like google earth.

%% Swath Width

Ec = 40075; % Earth Circumference in km

SwathWidth = Ec / Norbits;

% solving a right triangle for theta:

b = SwathWidth/2;
h = a - 6378;

ang = 2*(atan(b/h));

angleDeg = ang * (180/pi)