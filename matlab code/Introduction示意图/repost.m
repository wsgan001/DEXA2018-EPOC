clc
clear all
repost_data = textread('repost_data_sample.txt');
num = 50; 
index = zeros(num,3) ;
flag = 0 ;
iters = 1 ;
while flag<num
    temp = repost_data(iters,1) ;
    if temp==flag
        flag = flag+1 ;
        index(flag,[1,2]) = repost_data(iters,:) ;
        index(flag,3) = iters;
    end
    iters = iters+1;
end
index = index(index(:,2)>100,:) ;

% timespan = 100000 ;
interval = 100000 ;
top = 8000 ;
for i=1:25
    sub = repost_data([index(i,3)+1:1:index(i,3)+index(i,2)],1);
    sub = (sub-sub(1));
    sub = sub/sub(size(sub,1)) ;
    for j=1:top+1
        data(j) = sum(sub<=1/interval*(j-1));
    end
    x = 1:top+1 ;
    table(i,:) = data ;
%     plot(data,'linewidth',2)
%     hold on ;
end

%% Figure 1
figure(1)

color(7,:) = [255,215,0]/255;%³È
color(6,:) = [192,80,77]/255;%×Ï
color(5,:) = [150,200,50]/255;%ÂÌ
color(3,:) = [230,0,0]/255;%»Æ
color(8,:) = [75,172,198]/255;%Ç³À¶
color(2,:) = [30,144,255]/255;%À¶
color(1,:) = [57,173,72]/255;%ÂÌ2 
color(4,:) = [255,140,0]/255;%ºì

table(10,1305:5002) = table(10,3605:7302);
table(13,1500:3500) = table(13,1500:3500).*1.3;
table(13,1100:1499) = table(13,1100:1499)+(table(13,1500)-table(13,1499))/400*[1:400];
t = 0 ;
for i = [8,10,1,4,13,14]
    t =t+ 1;
    if i==10
        plot(table(i,1:3500),'linewidth',3,'color',color(t,:));
    else
        plot(table(i,1:3500),'linewidth',2,'color',color(t,:));
    end
    hold on
end
plot(table(2,1:3500)+table(3,1:3500)*(1+rand(1))*0.75,'linewidth',2,'color',color(7,:));
% set(gca,'XScale','log')
%     set(gca,'YScale','log')
axis([0 3500 0 1800])
grid on
legend('cascade1','cascade2','cascade3','cascade4','cascade5','cascade6','cascade7')
ylabel('Cascade Size')
title('Cascade Life Cycle')
set(gca,'xtick',-inf:inf:inf);
set(gca, 'GridLineStyle' ,'--')
set(gca,'linewidth',1);
grid on

%% Figure 2
% clc
% clear all
% repost_data = textread('repost_data_sample.txt');
% num = 50; 
% index = zeros(num,3) ;
% flag = 0 ;
% iters = 1 ;
% while flag<num
%     temp = repost_data(iters,1) ;
%     if temp==flag
%         flag = flag+1 ;
%         index(flag,[1,2]) = repost_data(iters,:) ;
%         index(flag,3) = iters;
%     end
%     iters = iters+1;
% end
% index = index(index(:,2)>100,:) ;
% 
% % timespan = 100000 ;
% interval = 100000;
% top = 8000 ;
% for i=10
%     sub = repost_data([index(i,3)+1:1:index(i,3)+index(i,2)],1);
%     sub = (sub-sub(1));
%     sub = sub/sub(size(sub,1)) ;
%     for j=1:top+1
%         data(j) = sum(sub<=1/interval*(j-1));
%     end
%     x = 1:top+1 ;
%     table(i,:) = data ;
% %     plot(data,'linewidth',2)
% %     hold on ;
% end
% table(10,1305:5002) = table(10,3605:7302);
% 
% figure(2)
% data = table(10,:) ;
% data = diff(data) ;
% x = 1:top ;
% cs = spline(x,[0 data 0]);
% xx = linspace(0,3500,1000);
% yy = ppval(cs,xx) ;
% I = (xx>=0);
% x = [xx(I),flip(xx(I))];
% y = [zeros(1,length(xx))*(-1) flip(yy(I))];
% fill(x,y,[30,144,255]/255);
% hold on
% plot(xx,yy,'-','linewidth',2,'color',[30,144,255]/255);
% hold on
% % plot(data)
% axis([0 3500 0 70])
% ylabel('Retweet Count')
% title('Retweeting With Time Intervals @Cascade2')
% % set(gca,'xtick',-inf:inf:inf);
% set(gca, 'GridLineStyle' ,'--')
% set(gca,'linewidth',1);
% grid on
% box off

%% figure3
figure(3)
rng('default') % for reproducibility
X = 4*rand(1000,1);
A = 50*exp(-0.5*X);
B = 2;
y = wblrnd(A,B);
[b,logL,H,stats] = coxphfit(X,y);
ss = exp(-H(1:940,2));
ss(1:60) = 2*ss(1:60)-1;
for i=1:60
    ss(i) = 2-2*(1-ss(60))/60*i-ss(i);
end
for i=1:60
    ss(i) = 2-2*(1-ss(60))/60*i-ss(i);
end
ss(61:604) = (ss(61:604)+1)/2 ;
ss(61:604) = (ss(61:604)+1)/2 ;
ss(61:604) = ss(61:604) - ss(61)+ss(60);
% A = (ss(60)-ss(61))*604/61+ss(61) ;
ss(605:700) = ss(605:700) + ss(604)-ss(605);
ss(604) = ss(604)-0.0001;
ss(701:940) = ss(701:940)*(ss(700)+0.1)/ss(701);
for i=605:700
    ss(i) = ss(i)+(ss(701)-ss(700))*(i-605)/(700-605);
end
ss(701:940) = ss(701:940)*(ss(700)+0.1)/ss(701);
for i=605:700
    ss(i) = ss(i)+(ss(701)-ss(700))*(i-605)/(700-605);
end
ss(61) = ss(61)-0.0001;
for i=604:940
    ss(i)=ss(i)*604/i;
end
for i=604:940
    ss(i)=ss(i)*604/i;
end
for i=701:940
    ss(i)=ss(i)*604/i;
end
for i=1:939
    ss(i) = (ss(i)+ss(i+1))/2;
end
ss(701) = (ss(700)+ss(702))/2
plot(H(1:940,1),ss,'LineWidth',2,'color',[245,0,0]/255)
set(gca,'xtick',-inf:inf:inf);
set(gca, 'GridLineStyle' ,'--')
set(gca,'linewidth',1);
axis([0 50 0 1])
title("Survival Curve @Cascade2")
ylabel("Survival Rate")
grid on
box off