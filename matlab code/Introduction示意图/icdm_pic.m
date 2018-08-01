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


a = [1:1:400];
table(10,1001:1400) = table(10,1001:1400)+0.001*a.^2;
table(10,1401:1450)= table(10,1401:1450)+[150:-149/49:1];
b=[0:1:350];
table(10,1000:1350) = table(10,1000:1350)+0.001*b.^2;
table(10,1351:1450)= table(10,1351:1450)+[150:-149/99:1];
for i = [8,10,1]
    t =t+ 1;
    plot(table(i,1:3500),'linewidth',3,'color',color(t,:));
    hold on
end

%plot(table(2,1:3500)+table(3,1:3500)*(1+rand(1))*0.75,'linewidth',2,'color',color(7,:));
% set(gca,'XScale','log')
%     set(gca,'YScale','log')
axis([0 3500 0 1800])
grid on
legend('cascade1','cascade2','cascade3')
ylabel('Cascade Size')
title('Cascade Life Cycle')
%set(gca,'xtick',-inf:inf:inf);
set(gca, 'GridLineStyle' ,'--')
set(gca,'linewidth',1);
grid on

