%% TOOL INPUTS
%10 dimensional vector
%clear all %helps to stop variables from getting mixed up when running the code

toolInput = [7000 0.01 50 0 0 0, 3 0 0 0];
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

%Plot Ground Track
figure(1)
plot(vecLon,vecLat,'.r');
title("Ground Track on Longitude Latitude Map");
xlabel("Longitude");
ylabel("Latitude");
grid on
hold on

%Set Image
img = imread('earthmap.jpg'); 
display = image(xlim,-[-90 90],img); 
uistack(display,'bottom')
axis([-180 180 -90 90])

grid on
hold off
%pause at the end of each section to give plots time to load
pause(1)
%% PROBLEM 2. Body-Fixed ECEF Trajectory w/ Scaled Earth

figure(2)
%Set Earth Image to Sphere
[x2,y2,z2] = sphere(40);
Earth = surf(x2 * radEarth, y2 * radEarth, z2 * radEarth);
img2 = flip(img,1);
set(Earth,'CData', img2,'FaceColor','texture');
hold on

%Plot 3D trajectory
plot3(rv(:,1), rv(:,2), rv(:,3),'-r');
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

figure(3)
%Set Earth Image to Sphere
[x3,y3,z3] = sphere(40);
Earth = surf(x3 * radEarth, y3 * radEarth, z3 * radEarth);
img3 = flip(img,1);
set(Earth,'CData', img3,'FaceColor','texture');
hold on

%Plot 3D trajectory
plot3(rBody(:,1), rBody(:,2), rBody(:,3),'-r');
axis equal
title('Inertial Frame Spacecraft Trajectory ')
xlabel('X km')
ylabel('Y km')
zlabel('Z km')
legend('Earth','Trajectory')
hold off

%pause at the end of each section to give plots time to load
pause(1)
%% Azimuth and Elevation vs. Time

%set variables
lon = 262.25;
lat = 30.3;
G = vecLon-lon;

%create elevation vector
for i=1:length(G)
    
Elevation(i) = atand((cosd(G(i))*cosd(lat)-0.1512)/(sqrt(1-((cosd(G(i))).^2)*((cosd(lat)).^2))));

end

%create azimuth vector
Azimuth = 180 + atand( tand(G)/sind(30.3) );

%
timeProb4 = 1:length(G);

figure(4)
plot(timeProb4,Azimuth,timeProb4,Elevation)
title('Azimuth and Elevation vs Time')
legend('Azimuth','Elevation')
xlabel('Time (min)')
ylabel('Degrees')
%pause at the end of each section to give plots time to load
pause(1)
%% Ground Track Animation

figure(5)

filename = 'AnimatedGroundTrackCaseA.gif';

%in steps of 5 to reduce file size and increase gif animation general speed
for j = 1:5:length(tspan)
 
plot(vecLon(1:j),vecLat(1:j),'.r');
hold on
xlim([-180 180]);
ylim([-90 90]);

img5 = imread('earthmap.jpg'); 
display5 = image([-180 180],-[-90 90],img5); 
uistack(display5,'bottom')
axis([-180 180 -90 90])
drawnow

title( 'Ground Track Animation ')
xlabel("Longitude");
ylabel("Latitude"); 

    frame = getframe(5);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if j==1
        imwrite(imind,cm,filename,'gif','Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end

end
hold off
pause(1)