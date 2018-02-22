clear
clc

filename = 'sandbed1.gif';
kk = 0;
for n = 0:1:500
    kk = kk + 1
    fig = figure(kk);
    s1 = num2str(n);
    s2 = '000.txt';
    s = [s1,s2];
    fid = fopen(s);
    a = textscan(fid, '%s%f%f%f%f%f');
    fclose(fid);    
    x0 = a(:,4);
    y0 = a(:,5);
    r0 = a(:,6);
    x = cell2mat(x0);
    y = cell2mat(y0);
    r = cell2mat(r0);
    alpha = 0:pi/30:2*pi;
    lengh = size(x);
	for i = 1:lengh
    	xx = x(i) + r(i)*cos(alpha);
   		yy = y(i) + r(i)*sin(alpha);
		if r(i)>0.49e-3
  			fill( xx, yy, '-k')
		elseif r(i)>0.41e-3
  			fill( xx, yy, '-b')
		elseif r(i)>0.39e-3
  			fill( xx, yy, '-r')
		elseif r(i)>0.31e-3
  			fill( xx, yy, '-g')
		elseif r(i)>0.29e-3
  			fill( xx, yy, '-y')
		elseif r(i)>0.21e-3
  			fill( xx, yy, '-c')
		elseif r(i)>0.19e-3
  			fill( xx, yy, '-w')
		elseif r(i)>0.11e-3
  			fill( xx, yy, '-p')
		else
  			fill( xx, yy, '-w')
		end
    	hold on
	end 
    	plot([0,0.045],[0,0],'linewidth',5,'color','k')
    	plot([0,0],[0,0.01],'linewidth',5,'color','k')
    	plot([0.045,0.045],[0,0.01],'linewidth',5,'color','k')
    	%plot([0,0.045],[0.045,0.045],'linewidth',5)
    	axis equal off;
    	drawnow
		%if kk<100
		%	axis( [0,0.045,0,0.045] )
		%else
			axis( [0,0.045,0,0.015] )
		%end 
     	set (gcf,'Position',[50,50,1266,568], 'color','w')
    	frame = getframe(kk);
    	im = frame2im(frame);
    	[A,map] = rgb2ind(im,256); 
    	if kk == 1;
        	imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.05);
    	else
        	imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.05);
    	end
    close(figure(kk));
end 
