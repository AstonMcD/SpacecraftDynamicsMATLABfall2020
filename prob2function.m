function drdt = prob2function(t,r,G,M1,M2,M3)

r12 = sqrt((r(1)-r(7))^2 + (r(2)-r(8))^2 + (r(3)-r(9))^2);
r13 = sqrt((r(1)-r(13))^2 + (r(2)-r(14))^2 + (r(3)-r(15))^2);
r23 = sqrt((r(13)-r(7))^2 + (r(14)-r(8))^2 + (r(15)-r(9))^2);

%VELOCITY 1%
drdt(1) = r(4);     %x
drdt(2) = r(5);     %y
drdt(3) = r(6);     %z
%ACCELERATION 1%
drdt(4) = 1*G*M2*(r(7)-r(1))/r12^3 + 1*G*M3*(r(13)-r(1))/r13^3; %x
drdt(5) = 1*G*M2*(r(8)-r(2))/r12^3 + 1*G*M3*(r(14)-r(2))/r13^3; %y
drdt(6) = 1*G*M2*(r(9)-r(3))/r12^3 + 1*G*M3*(r(15)-r(3))/r13^3; %z

%VELOCITY 2%
drdt(7) = r(10);      %x
drdt(8) = r(11);      %y
drdt(9) = r(12);      %z
%ACCELERATION 2%
drdt(10) = 1*G*M1*(r(1)-r(7))/r12^3 + 1*G*M3*(r(13)-r(7))/r23^3; %x
drdt(11) = 1*G*M1*(r(2)-r(8))/r12^3 + 1*G*M3*(r(14)-r(8))/r23^3; %y
drdt(12) = 1*G*M1*(r(3)-r(9))/r12^3 + 1*G*M3*(r(15)-r(9))/r23^3; %z

%VELOCITY 3%
drdt(13) = r(16);      %x
drdt(14) = r(17);      %y
drdt(15) = r(18);      %z
%ACCELERATION 3%
drdt(16) = 1*G*M1*(r(1)-r(13))/r13^3 + 1*G*M2*(r(7)-r(13))/r23^3; %x
drdt(17) = 1*G*M1*(r(2)-r(14))/r13^3 + 1*G*M2*(r(8)-r(14))/r23^3; %y
drdt(18) = 1*G*M1*(r(3)-r(15))/r13^3 + 1*G*M2*(r(9)-r(15))/r23^3; %z

drdt = drdt';
end