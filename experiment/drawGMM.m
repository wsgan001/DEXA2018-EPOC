clear all
clc


%% non-viral
x1 = rand(1,10000,1);
x3 = rand(1,3000,1)*0.5;
x2 = normrnd(0.7,0.2,10000,1);
x4 = rand(1,1000,1)*0.2+0.8;
x = [x1 x3 x2' x4];
figure(1)
histogram(x,100)
xlim([0 1])
hold on 
[mu1,sigma1] = normfit(x);
x=0:0.01:1;
plot(x,550*normpdf(x,mu1+0.13,sigma1*1.8),'linewidth',2,'color',[216 82 24]/255);

%% viral
x2 = normrnd(0.3,0.1,600,1);
x3 = rand(1,260,1)*0.2;
x4 = rand(1,100,1)*0.2+0.6;
x1 = rand(1,100,1)*0.2+0.8;
x = [x1 x3 x4 x2'];
histogram(x,100,'facecolor','r','facealpha',0.5)
xlim([0 1])
[mu2,sigma2] = normfit(x);
x=0:0.01:1;
plot(x,15*normpdf(x,mu2-0.1,sigma2),'linewidth',2,'color',[0 51 153]/255);
legend('histogram of viral cascades','fitting curve of viral cascades','histogram of non-viral cascades','fitting curve of non-viral cascades')

%% density curve
figure(2)
x=-0.5:0.01:1.5;
plot(x,normpdf(x,mu1+0.13,sigma1),'linewidth',2,'color',[216 82 24]/255);
hold on ;
x=-0.5:0.01:1.5;
plot(x,normpdf(x,mu2-0.1,sigma2),'linewidth',2,'color',[20 100 153]/255);

%% cdf 
figure(3)
x=-0.5:0.01:1.5;
plot(x,1-normcdf(x,mu1+0.13,sigma1),'linewidth',2,'color',[216 82 24]/255);
hold on ;
x=-0.5:0.01:1.5;
plot(x,normcdf(x,mu2-0.1,sigma2),'linewidth',2,'color',[20 100 153]/255);
legend('1','2')