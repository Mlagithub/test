clear
clc
a = load('0.txt');
xi = a(:,1);
x = a(:,2);
y = a(:,3);
r = a(:,4);
Br = max(x+r);
Bt = max(y)*1.2;
alpha = 0:pi/30:2*pi;
lengh1 = size(x,1);
lengh2 = size(alpha,2);
xx = zeros(lengh1,lengh2);
yy = zeros(lengh1,lengh2);
for i = 1:lengh1
    xx(i,:) = x(i) + r(i)*cos(alpha);
    yy(i,:) = y(i) + r(i)*sin(alpha);
end
plot(xx', yy','k')
hold on
plot([0,Br],[0,0],'linewidth',5,'color','k')
plot([0,0],[0,Bt],'linewidth',5,'color','k')
plot([Br,Br],[0,Bt],'linewidth',5,'color','k')
axis equal
hold on 
axis off
axis( [0,Br,0,Bt] )

%text(x(i), y(i), num2str(xi(i)));

% v = 0.17;
% G = 31.2e9;
% deltat = pi*(min(r)*sqrt(2650/G)/(0.163*v+0.877))
