clear
clc

Rho = 2650;
H_sandbed = 6.0e-3;
Length = 0.02;
CellE = 5e-3;
X3 = [0.784,1,1.038,1.015,0.969,0.918,0.867,0.818];
Rh = [10,20,30,40,50,60,70,80];

a = zeros(2000,6);
Q = zeros( int16((0.5-H_sandbed)/CellE),1 );
E = zeros( size(Q,1),1 );
H = zeros( size(Q,1),1 );
Field = zeros( 8,size(Q,1));




start = 1140;% 140 240--0.6£¬340--3£¬440--2.3£¬540--0.5£¬640--1.5£¬740--0.4,840--1£¬940--1.5
            % 1040--1.6 1140--2.7  1240--1.5  1340--1.7 
nfile = 1;
kk = 0;
for n = start+1:1:start+nfile*100
    kk = kk + 1
    s1 = num2str(kk+start);
    s2 = '.txt';
    s = [s1,s2];
    a = load(s);
    b = find(a(:,6)~=0);
    for i=1:size(b,1)
        j = int16( (a(b(i),3)-H_sandbed)/CellE )+1;
        Q(j) = Q(j) + a(b(i),6);
    end
end
Q = Q*1e-2/nfile;


for i = 1:size(Q,1)
    H(i) = CellE*i;
    for j = 1:size(Q,1)
        if( j~=i )
            distance = CellE*(j-i);
            E(i) = E(i) + 9e9*Q(j)/distance^2;
        end
    end
end


for i = 1:8
    Field(i,:) = E.*X3(i);
end


% % Electric Field with height
% figure(1)
% % semilogy(H,E/1e3,'linewidth',5,'color','b')
% % subplot( 2,2,4 )
% plot(H,E/1e3,'linewidth',5,'color','b')
% legend( {'Rh=20% Temperature=20^{o}C'},'FontSize',16,'FontWeight','bold' )
% xlabel('Height (m)','FontSize',16,'FontWeight','bold','Color','k')
% ylabel('Electric Field Intensity(Kv/m)','FontSize',16,'FontWeight','bold','Color','k')
% set( gca,'Fontsize',14 )
% text( 0.25*max(H(:,1)),0.5*max(E(:,1))/1e3,'Time=10(s),Height=0.05(m)','FontSize',16,'FontWeight','bold' )
% hold on

% Electric at some point 
figure(2)
plot( Rh,max(Field,[],2)/0.674e3,'color','b','linewidth',4,'marker','+','markeredgecolor','r','markersize',10 )
% legend( {'U_{*}=0.4(m/s)'},'FontSize',16,'FontWeight','bold' )
hold on
plot( Rh,max(Field,[],2)/1.0e3,'color','r','linewidth',4,'marker','+','markeredgecolor','r','markersize',10 )
% legend( {'U_{*}=0.5(m/s)'},'FontSize',16,'FontWeight','bold' )
hold on
plot( Rh,max(Field,[],2)/1.365e3,'color','g','linewidth',4,'marker','+','markeredgecolor','r','markersize',10 )
legend( {'U_{*}=0.6(m/s)','U_{*}=0.5(m/s)','U_{*}=0.4(m/s)'},'FontSize',16,'FontWeight','bold' )
hold on
xlabel('Rh(%)','FontSize',16,'FontWeight','bold','Color','k')
ylabel('Max Electric Field Intensity(kv/m)','FontSize',16,'FontWeight','bold','Color','k')
% legend( {'U_{*}=0.5(m/s)'},'FontSize',16,'FontWeight','bold' )
text( 30,1.5,'Height=0.05(m)','FontSize',16,'FontWeight','bold' )
set( gca,'Fontsize',14 )
hold on
xlim([10,80]);
ylim([1.0,5]);
