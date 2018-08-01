clc
clear all
rng('default') % for reproducibility
X = 4*rand(100,1);
A = 50*exp(-0.5*X);
B = 1.9;
y = wblrnd(A,B);
% figure(1)
% [b,logL,H,stats] = coxphfit(X,y);
% 
% % stairs(H(:,1),1-(1-exp(-H(:,2)))/1.5,'LineWidth',2)
% % hold on ;
% % stairs(H(:,1),1-(1-exp(-H(:,2)))/1.2,'LineWidth',2)
% % stairs(H(:,1),1-(1-exp(-H(:,2)))/1.3,'LineWidth',2)
% % stairs(H(:,1),1-(1-exp(-H(:,2)))/1.8,'LineWidth',2)
% % stairs(H(:,1),1-(1-exp(-H(:,2)))/1.1,'LineWidth',2)
% % stairs(H(:,1),1-(1-exp(-H(:,2)))/1.7,'LineWidth',2)
xx = linspace(0,100);
% m = 1-wblcdf(xx,50*exp(-0.5*mean(X));
% m = 1-(1-m)/2 ;
yy = 1-wblcdf(xx,50*exp(-0.5*mean(X)),B)./1.4;
% line(xx,yy,'linestyle','--','color','r','LineWidth',2)
% hold on 
% 
% %% fig1
% % h1=ezplot('x-30-5*exp(-500*(y-0.40)^2)',[0 50 0.25 0.55]);
% % h2=ezplot('x-30+y*0.00001',[29 31 0 0.8]);
% % set(h2,'linestyle','--','color','black');
% y = 0.25:0.004:0.55;
% x = 30+3.5*exp(-600*(y-0.40).^2);
% plot(x,y,'color',[50,100,200]/255)
% [~,N]=size(y);
% for i =1:N-1
%     fill([x(i),x(i+1),30,30],[y(i),y(i+1),y(i+1),y(i)],[100,180,255]/255,'edgealpha',0,'FaceAlpha',0.4);
%     hold on
% end
% 
% y = 0.15:0.004:0.4;
% x = 30+3.5*exp(-1000*(y-0.27).^2);
% plot(x,y,'color',[50,100,200]/255)
% [~,N]=size(y);
% for i =1:N-1
%     fill([x(i),x(i+1),30,30],[y(i),y(i+1),y(i+1),y(i)],[230,0,0]/255,'edgealpha',0,'FaceAlpha',0.4);
%     hold on
% end
% y = 0:0.001:1.1;
% x = y.^0*30;
% plot(x,y,'--','color',[50,100,200]/255)
% 
% %% fig2
% % h1=ezplot('x-30-5*exp(-500*(y-0.40)^2)',[0 50 0.25 0.55]);
% % h2=ezplot('x-30+y*0.00001',[29 31 0 0.8]);
% % set(h2,'linestyle','--','color','black');
% y = 0.2:0.004:0.5;
% x = 40+3.7*exp(-400*(y-0.38).^2);
% plot(x,y,'color',[50,100,200]/255)
% [~,N]=size(y);
% for i =1:N-1
%     fill([x(i),x(i+1),40,40],[y(i),y(i+1),y(i+1),y(i)],[100,180,255]/255,'edgealpha',0,'FaceAlpha',0.4);
%     hold on
% end
% 
% y = 0.1:0.004:0.35;
% x = 40+3.7*exp(-700*(y-0.23).^2);
% plot(x,y,'color',[50,100,200]/255)
% [~,N]=size(y);
% for i =1:N-1
%     fill([x(i),x(i+1),40,40],[y(i),y(i+1),y(i+1),y(i)],[230,0,0]/255,'edgealpha',0,'FaceAlpha',0.4);
%     hold on
% end
% y = 0:0.001:1.1;
% x = y.^0*40;
% plot(x,y,'--','color',[50,100,200]/255)
% 
% %% fig3
% % h1=ezplot('x-30-5*exp(-500*(y-0.40)^2)',[0 50 0.25 0.55]);
% % h2=ezplot('x-30+y*0.00001',[29 31 0 0.8]);
% % set(h2,'linestyle','--','color','black');
% y = 0.8:0.004:1;
% x = 7.5+2.5*exp(-1500*(y-0.9).^2);
% plot(x,y,'color',[50,100,200]/255)
% [~,N]=size(y);
% for i =1:N-1
%     fill([x(i),x(i+1),7.5,7.5],[y(i),y(i+1),y(i+1),y(i)],[100,180,255]/255,'edgealpha',0,'FaceAlpha',0.4);
%     hold on
% end
% 
% y = 0.75:0.004:0.9;
% x = 7.5+2.5*exp(-2800*(y-0.825).^2);
% plot(x,y,'color',[50,100,200]/255)
% [~,N]=size(y);
% for i =1:N-1
%     fill([x(i),x(i+1),7.5,7.5],[y(i),y(i+1),y(i+1),y(i)],[230,0,0]/255,'edgealpha',0,'FaceAlpha',0.4);
%     hold on
% end
% y = 0:0.001:1.1;
% x = y.^0*7.5;
% plot(x,y,'--','color',[50,100,200]/255)
% 
% xlim([0,50])


%% hazadous function
figure(3)
y1 = [0 -diff(yy)]./yy;
y1(1:18) = (y1(1:18)+y1(18))/2 ;
line(xx(1:45),y1(1:45),'color',[1 0 0],'linestyle','--','LineWidth',2);
% line(xx(1:45),y1(1:45)+0.001,'color',[0 0 0],'linestyle','--','LineWidth',2);
% line(xx(1:45),y1(1:45)-0.001,'color',[0 0 0],'linestyle','--','LineWidth',2);
hold on
[~,N]=size(y1);
for i =1:44
    fill([xx(i),xx(i+1),xx(i+1),xx(i)],[y1(i)-0.003,y1(i+1)-0.003,y1(i+1)+0.003,y1(i)+0.003],[0 0 0],'edgealpha',0,'FaceAlpha',0.05);
    hold on
end
xlim([0,50])
ylim([0,0.06])
legend('baseline hazard rate')
% X = 4*rand(100,1);
% A = 50*exp(-0.5*X);
% B = 1.4;
% y = wblrnd(A,B);
% xx = linspace(0,100);
% yy = 1-wblcdf(xx,50*exp(-0.5*mean(X)),B)./1.4;
% y2 = [0 -diff(yy)]./yy;
% y2(1:15) = (y2(1:15)+2*y2(15))/3;
% y2(1:20) = (0.05+y2(1:20))/2;
% for i =20:51
%     y2(i) = y2(i)+0.0001*(i-35)^2;
% end
% t = y2(40);
% for i =20:51
%     y2(i) = y2(i)-t;
% end
% y2(19) = y2(19)-0.0003;
% y2(21) = y2(21)-0.00005;
% y2(20)=0.0425;
% y2(21:51) = (5*y2(21:51)+y2(21))/6;
xx = 30:0.01:51;
y2 = 0.03*exp(-0.0022*(xx-30).^2);
xx2 = 0:0.01:30 ;
y3 = 0.03*exp(-0.001*(xx2-30).^2);
xx = [xx2 xx];
y2 = [y3 y2];
line(xx(1:4450),y2(1:4450),'color',[250,100,100]/255,'LineWidth',2);

B = 3;
y = wblrnd(A,B);
xx = linspace(0,100);
yy = 1-wblcdf(xx,50*exp(-0.5*mean(X)),B)./1.4;
y2 = [0 -diff(yy)]./yy;
y2(1:20) = (y2(1:20)*10+y2(20))/11;
y2(20:51) = y2(20:51) +0.001;
line(xx(1:45),y2(1:45)/6,'color',[40,200,240]/255,'LineWidth',2);

% ylim([0,6])

% y = normrnd(0.65,0.1,100);
% [~,N]=size(y);
% figure(2)
% m = histogram(y,100) ;
% figure(1)
% x = m.BinEdges;
% y = m.Values/200 ;
% plot(x(1:end-1),y);
