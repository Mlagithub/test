clear
clc

Br = 0.02;
Bt = 0.07;
%Br = 6e2;
%Bt = 6e2;
filename = '01.gif';
start = 0;
kk = start;
nfile = input('Input the number of files: ');
for n = 1:1:nfile
    kk = kk + 1
    fig = figure(kk);
    s1 = num2str(kk);
    s2 = '.txt';
    s = [s1,s2];
    a = load(s);
    xi = a(:,1);
    x = a(:,2);
    y = a(:,3);
    r = a(:,4);
    alpha = 0:pi/30:2*pi;
    lengh1 = size(x,1);
    lengh2 = size(alpha,2);
    xx = zeros(lengh1,lengh2);
    yy = zeros(lengh1,lengh2);
    ss = zeros(lengh1,3);
    for i = 1:lengh1
        xx(i,:) = x(i) + r(i)*cos(alpha);
        yy(i,:) = y(i) + r(i)*sin(alpha);
%         if r(i)==1
%             ss(i,:) = [1,0,0];
%         elseif r(i)==2
%             ss(i,:) = [0,1,0];
%         elseif r(i)==3
%             ss(i,:) = [0,0,1];
%         else
%             ss(i,:) = [0.1,0.8,0.2];
%         end 
%         fill(xx(i,:),yy(i,:),ss(i,:))
%         hold on
    end
    plot(xx', yy')
    hold on
	plot([0,Br],[0,0],'linewidth',5,'color','k');
	plot([0,0],[0,Bt],'linewidth',5,'color','k');
	plot([Br,Br],[0,Bt],'linewidth',5,'color','k')
	set(gcf,'visible','off')
	axis equal ;
    %axis off
	drawnow;
	axis( [0,Br,0,Bt] )
	set (gcf,'Position',[1,41,1366,650], 'color','w')
	frame = getframe(kk);
	im = frame2im(frame);
	[A,map] = rgb2ind(im,256); 
	if kk == start+1;
		imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.5);
	else
		imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.5);
	end
    close(figure(kk));
end 


