function y = FiveCompFit
T = readmatrix('FiveCompData.xlsx');
ix = [0, 2, 5, 15, 30, 60, 120, 240, 480, 1440, 2880];

lowerBond = [0, 0, 0, 0, 0, 0, 0, 0, 0];
upperBond = [1, 1, 1, 1, 1, 1, 1, 1, 1];

OLS = @(b)FiveCompOLS(b);

B = particleswarm(OLS, 9, lowerBond, upperBond);

odefun = @(t,y,b) [-(b(1)+b(3)+b(5)+b(9))*y(1)+b(2)*y(2)+b(4)*y(3)+b(6)*y(4);
-(b(2))*y(2)+b(1)*y(1);
-(b(4))*y(3)+b(3)*y(1);
-(b(6)+b(7))*y(4)+b(5)*y(1)+b(8)*y(5);
-b(8)*y(5)+b(7)*y(4)];


y0 = T(1, 5:9);
[~, y] = ode45(@(t,Y)odefun(t,Y,B), ix, y0); 
end
