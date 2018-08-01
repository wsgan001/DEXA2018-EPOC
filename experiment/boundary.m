cumhaz = xlsread('cumhaz');
% plot(exp(-cumhaz))
b = [4.852*10^-2  8.195*10^-6 -5.629*10^-3] ;
ht = diff(cumhaz);
X1 = [] ;% cascade_num
X2 = [] ;% follower
X3 = [] ;% original timestamp
for i =1:cas_num
    i
    for j=1:round(max(cascade_list(i,:))/interval)
        X1(i,j) = length(find(cascade_list(i,:)<=j*interval & cascade_list(i,:)>=0)) ;
        X2(i,j) = sum(follower_list(i,(find(cascade_list(i,:)<=j*interval & cascade_list(i,:)>=0 )))) ;
        X3(i,j) = index(i,1);
    end
end

interval = 60;
for i=1:1000
    for j = 1:168
        if X1(i,j)==0
            break
        end
    end
    if j~=1
        X1(i,j:168)=X1(i,j-1)*ones(1,168-j+1);
    end
end
for i=1:1000
    for j = 1:168
        if X2(i,j)==0
            break
        end
    end
    if j~=1
        X2(i,j:168)=X2(i,j-1)*ones(1,168-j+1);
    end
end
for i=1:1000
    for j = 1:168
        if X3(i,j)==0
            break
        end
    end
    if j~=1
        X3(i,j:168)=X3(i,j-1)*ones(1,168-j+1);
    end
end

cumht(1:1000,1:4529)=0;
for i=1:1000
    for j=2:4529
        cumht(i,j) = cumht(i,j-1)+ht(j)*exp(b(1)*X1(i,floor(j/interval)+1)+b(2)*X2(i,floor(j/interval)+1)+b(3)*X3(i,floor(j/interval)+1));
    end
end        

% for i=1:1000
%     if index(i,4)==1
%         plot(exp(-cumht(i,:)),'b');
%         hold on ;
%     else
%         plot(exp(-cumht(i,:)),'r');
%         hold on;
%     end
% end
% cumht = exp(-cumht);
% for i=1:4529
%     tempV = cumht(find(index(:,4)==1 ),i);
%     tempN = cumht(find(index(:,4)==0 ),i);
%     [mu1,sigma1] = normfit(tempV);
%     [mu2,sigma2] = normfit(tempN);
%     S(i) = (mu1*max(sigma2,0.00001)+mu2*max(sigma1,0.00001))/(max(sigma1,0.00001)+max(sigma2,0.00001));
% end
% hold on ;
% judge = zeros(1,200);
% for i  =801:1000
%     for j =2:2
%         if cumht(i,j)<S(j)
%             judge(i-800) = 1;
%             break;
%         end
%     end
% end

length(find(index(801:1000,4)==1 & judge'==1))/length(find((judge'==1)));
% % tempV
% plot(0:0.01:1,normpdf(0:0.01:1,mu1,sigma1))
% hold on ;
% figure(3)
% histogram(tempV,100)
% % tempN
% plot(0:0.01:1,normpdf(0:0.01:1,mu2,sigma2))
% hold on ;
% figure(3)
% histogram(tempN,100)

% plot(cumht(find(cumht(:,100)>0.1),:);

%% a 
ok = zeros(1000,4529);
index2 = zeros(1000,1);
t = 1;
for i=1:cas_num
    temp = cumht(i,:);
    flag = 0;
    for j=1:4528
        if temp(j)-temp(j+1)>0.01
            flag = 1;
            break;
        end
    end
    if flag == 0
        ok(t,:) = temp;
        index2(t) = 1-2*index(i,4) ;
        t = t+1;
    end
end
% [linspace(0,0.2,2264),linspace(0.2,0,2265)]'+
figure(2)
[~,index_index] = sort(ok(1:120,1000));
plot(ok(170,:)','color',[100,180,220]/255);
hold on;
plot(ok(141,:)','color',[255,0,0]/255);
plot((ok(170:180,:)'),'color',[100,180,220]/255);
plot((ok(140:145,:)'),'color',[255,0,0]/255);
ok = ok(index_index,:);
plot((ok(25:120,:)'),'color',[100,180,220]/255);
plot(ok(1:25,:)','color',[255,0,0]/255);
plot(ok([12,13],:)','color',[100,180,220]/255);
legend('non-viral','viral')
xlim([1 1800])
xlabel('Minute')
ylabel('Surviral rate')
title('Cox-PH time-dependent');

