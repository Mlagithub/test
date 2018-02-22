clear
clc


filename = 'sandbedfl00.gif';
	h = figure(1);
    fid = fopen('fl.txt');
    a = textscan(fid, '%s%f%f%f');
    fclose(fid);    
    x0 = a(:,2);
    y0 = a(:,3);
    r0 = a(:,4);
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
    	plot([0,0.045],[0,0],'linewidth',5,'color','k')
    	plot([0,0],[0,0.01],'linewidth',5,'color','k')
    	plot([0.045,0.045],[0,0.01],'linewidth',5,'color','k')
    	%plot([0,0.045],[0.045,0.045],'linewidth',5)
    	axis equal off;
    	drawnow
		axis( [0,0.045,0,0.015] )
    	set (gcf,'Position',[1,41,1366,650], 'color','w')
    	frame = getframe(h);
    	im = frame2im(frame);
    	[A,map] = rgb2ind(im,256); 
    	if i == 1;
        	imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.05);
    	else
        	imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.05);
    	end
	end 
    close(figure(1));
