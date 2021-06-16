function LSD = FiveCompOLS(b)
tspan = [0, 2, 60, 120, 240, 480, 1440, 2880];
tblood = [0, 5, 15, 30, 60];
T = readmatrix('FiveCompData.xlsx');
yxblood = T(1:5,1); yxkidney = T(1:8,2); yxliver = T(1:8,3); yxhead = T(1:8,4);
ix = [0, 2, 5, 15, 30, 60, 120, 240, 480, 1440, 2880];
iyxblood = interp1(tblood, yxblood, ix)';
iyxkidney = interp1(tspan, yxkidney, ix)';
iyxliver = interp1(tspan, yxliver, ix)';
iyxhead = interp1(tspan, yxhead, ix)';


odefun = @(t,y,b) [-(b(1)+b(3)+b(5)+b(9))*y(1)+b(2)*y(2)+b(4)*y(3)+b(6)*y(4);
-(b(2))*y(2)+b(1)*y(1);
-(b(4))*y(3)+b(3)*y(1);
-(b(6)+b(7))*y(4)+b(5)*y(1)+b(8)*y(5);
-b(8)*y(5)+b(7)*y(4)];

y0 = T(1, 5:9);

[~, y] = ode45(@(t,Y)odefun(t,Y,b), ix, y0);

LSD = sum(((y(:,2) - iyxkidney)/max(iyxkidney)).^2)+sum(((y(:,3) - iyxliver)/max(iyxliver)).^2) ...
    +sum(((y(:,4) - iyxhead)/max(iyxhead)).^2)+sum(((y(:,1) - iyxblood)/max(iyxblood)).^2);
end