clear
clc

r = 90e-6;
miv = 0.3;
E = 100e6;
G = E/2/(1+miv);
deltat = pi*r*sqrt(2650/G)/(0.163*miv+0.887);


%a = load('SandBed');
a = load('400.txt');

figure(1)
xi = a(:,1);
x = a(:,2);
y = a(:,3);
r = a(:,4);
Br = 0.02;
Bt = 0.01;
alpha = 0:pi/30:2*pi;
lengh1 = size(x,1);
lengh2 = size(alpha,2);
xx = zeros(lengh1,lengh2);
yy = zeros(lengh1,lengh2);
ss = zeros(lengh1,3);
for i = 1:lengh1
    xx(i,:) = x(i) + r(i)*cos(alpha);
    yy(i,:) = y(i) + r(i)*sin(alpha);
%     if r(i)==1
%         ss(i,:) = [1,0,0];
%     elseif r(i)==2
%         ss(i,:) = [0,1,0];
%     elseif r(i)==3
%         ss(i,:) = [0,0,1];
%     else
%         ss(i,:) = [0.1,0.8,0.2];
%     end 
%     fill(xx(i,:),yy(i,:),ss(i,:))
%     hold on
end
plot(xx', yy')
%text(x(:,1),y(:,1),num2str(xi(:,1)))
hold on
plot([0,Br],[0,0],'linewidth',5,'color','k')
plot([0,0],[0,Bt],'linewidth',5,'color','k')
plot([Br,Br],[0,Bt],'linewidth',5,'color','k')
axis equal
hold on 
%axis off
%axis( [0.5*Br,0.55*Br,0,Bt] )
%axis( [0,Br,0,Bt] )
%text(x(i), y(i), num2str(xi(i)));


% a = load('log.txt');
% x = a(:,1);
% y = a(:,2);
% figure(2)
% plot(x,y,'o');
