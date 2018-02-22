clear
clc

Rho = 2650;
nfile = 2000;

a = zeros(2000,6);
Mass = zeros(nfile,1);


start = 300;
kk = 0;
for n = start+1:1:nfile
    kk = kk + 1
    s1 = num2str(kk+start);
    s2 = '.txt';
    s = [s1,s2];
    a = load(s);
    b = find(a(:,6)~=0);
    for i=1:size(b,1)
        Mass(kk) = Mass(kk) + Rho*4*pi*a(b(i),4)^3/3;
    end
end

% load Mass.mat
T = [1e-2:1e-2:20];
figure(1)
plot(T,Mass/1e-2,'color','b','linewidth',2)
xlabel('Time (s)','FontSize',18,'FontWeight','bold','Color','k')
ylabel('Total mass of sand (kg)','FontSize',18,'FontWeight','bold','Color','k')
legend( {'Total mass of sand'},'FontSize',16,'FontWeight','bold' )
set( gca,'Fontsize',14)



