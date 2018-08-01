%% ############################################# load data #######################################################
clear all 
clc
% data = textread('data.txt');
% index = textread('index.txt');
data = textread('data_small.txt');
index = textread('index_small.txt');
index = index(1:1000,:);
data = data(1:index(1000,3),:);
% index = index2(1:100,:) ;
[cas_num,~] = size(index);
% data = data2(1:index(cas_num,3),:);
[tweet_num,~]  =size(data);

Train = 0.8 ;
%% ###################################################### list cascade ############################################
% get list
cascade_list = ones(cas_num,max(index(:,3)-index(:,2))+1)*(-1) ;
follower_list = ones(cas_num,max(index(:,3)-index(:,2))+1)*(-1) ;
t = 0 ;
cascade_index = 0 ;
for i=1:tweet_num
    if i==1
        cascade_index = cascade_index + 1 ;
        t = 1;
        cascade_list(cascade_index,t) = data(i,1) ;
        follower_list(cascade_index,t) = data(i,2) ;
    else
        if data(i,1)==0 && data(i-1,1)~=0
            cascade_index = cascade_index + 1 ;
            t = 1;
            cascade_list(cascade_index,t) = data(i,1) ;
            follower_list(cascade_index,t) = data(i,2) ;
        else
            t = t+1 ;
            cascade_list(cascade_index,t) = data(i,1) ;
            follower_list(cascade_index,t) = data(i,2) ;
         end
    end
end
cascade_list = cascade_list/60 ; % 分钟
% follower_list = cumsum(follower_list')';

%% ############################################random sort############################################
shuffle = randperm(cas_num);
cascade_list = cascade_list(shuffle,:);
follower_list = follower_list(shuffle,:);
index = index(shuffle,:);
%% ############################################ get threshold ############################################
sum_rank = sort(index(1:cas_num*Train,3)-index(1:cas_num*Train,2)) ;
threshold = sum_rank(round(cas_num*Train*0.95));
index(:,4) = (index(:,3)-index(:,2))>threshold ;

for i =1:cas_num
    temp = cascade_list(i,threshold);
    if temp>0
        burst_time(i) = temp;
    else
        burst_time(i) = -1;
    end
end

%% ############################################ get feature: total_sum, total_folower, orignial time, Ylabel ############################################
X = zeros(cas_num,4);
time_window = 9;
hpara = [0.503840700176885;1.34891941210983;1.61646027474011;0.0902052036239316;0.307628889573712;0.374075369154548];
npara = [1.16034926111027;24.5113780989613;0.892817833917356;0.354705514409481;9.59764524619549;7.09710678252588] ;
for epoch=1:round(max(max(cascade_list))/60)
    fprintf('epoch %d start\n',epoch);
    observe_time = 60*epoch; %一小时
    for i =1:cas_num
        X(i,1) = length(find(cascade_list(i,:)<observe_time & cascade_list(i,:)>=0));
        X(i,2) = sum(follower_list(find(cascade_list(i,:)<observe_time & cascade_list(i,:)>=0)));
        X(i,3) = index(i,1);
        X(i,4) = index(i,3)-index(i,2);
    end
    for i=1:cas_num
        temp = cascade_list(i,find(cascade_list(i,:)<observe_time & cascade_list(i,:)>=0));
        testX(i,1) = length(temp);
        for j = 1:time_window
            temp2(j) = length(temp(find(temp<=j*temp(end)/time_window)));
        end
        for j =1:time_window-1
            if temp2(j)~=0
                temp3(j) = (temp2(j+1)-temp2(j))/temp2(j);
            else
                temp3(j) = 0;
            end
        end
        testX(i,2) = max(sum(temp3)/length(temp3),0.01);
        testX(i,3) = sqrt(sum((temp2-sum(temp2)/length(temp2)).^2)/time_window);
        testX(i,4) = (index(i,3)-index(i,2))>=threshold;
    end
    testX(:,1) = testX(:,1)/(max(testX(:,1))+1) ;

    % SVR
    baseline1 = fitrsvm(X(1:cas_num*Train,1:3),X(1:cas_num*Train,4),'KernelFunction','gaussian','Standardize',true);
    predict_count1 = predict(baseline1,X(cas_num*Train+1:end,1:3));
    predictlabel1(1:round(cas_num*(1-Train)),epoch) = (predict_count1>=threshold);
    % Linear
    baseline2 = fitlm(X(1:cas_num*Train,1:3),X(1:cas_num*Train,4));
    predict_count2 = predict(baseline2,X(cas_num*Train+1:end,1:3));
    predictlabel2(1:round(cas_num*(1-Train)),epoch) = (predict_count2>=threshold);
    % Prewhether
    if epoch ==-1
        options = optimoptions('fmincon','Display','off');
        [hpara,~] = fmincon(@(para) likelihood_b1_1(testX(find(testX(1:cas_num*Train,4)==1),1:3),para),rand(6,1),[],[],[],[],[0,0,0,0,-200,0],[200,200,200,200,200,200], [] ,options);
        [npara,~] = fmincon(@(para) likelihood_b1_1(testX(find(testX(1:cas_num*Train,4)==0),1:3),para),rand(6,1),[],[],[],[],[0,0,0,0,-200,0],[200,200,200,200,200,200], [],options);
        fprintf('PreWhether finished!\n');
    end
    po_c(:,1) = length(find(testX(1:cas_num*Train,4)==1))/(cas_num*Train) .* betapdf(testX(cas_num*Train+1:end,1),hpara(1),hpara(2)) .* gampdf(testX(cas_num*Train+1:end,2),hpara(3),hpara(4)) .* normpdf(testX(cas_num*Train+1:end,3),hpara(5),hpara(6));
    po_c(:,2) = length(find(testX(1:cas_num*Train,4)==0))/(cas_num*Train) .* betapdf(testX(cas_num*Train+1:end,1),npara(1),npara(2)) .* gampdf(testX(cas_num*Train+1:end,2),npara(3),npara(4)) .* normpdf(testX(cas_num*Train+1:end,3),npara(5),npara(6));
    pco = [po_c(:,1) ./ (po_c(:,1) + po_c(:,2)) , po_c(:,2) ./ (po_c(:,1) + po_c(:,2))];
    predictlabel3(1:round(cas_num*(1-Train)),epoch) = (pco(:,1)>=pco(:,2));
end   
    

svr_time = ones(round(cas_num*(1-Train)),1)*(-1);
for i=1:round(cas_num*(1-Train))
    for j =1:round(max(max(cascade_list))/60)
        if predictlabel1(i,j)==index(i+cas_num*Train,4)
            svr_time(i,1)=j*60;
            break;
        end
    end
end
for i=1:round(cas_num*(1-Train))
    if burst_time(i+cas_num*Train)~=-1
        if svr_time(i,1)==-1
            svr_time(i,1)=0 ;
        else
            svr_time(i,1) = max((burst_time(i+cas_num*Train)-svr_time(i,1)),0);
        end
    else
        svr_time(i,1)=-1;
    end
end

linear_time = ones(round(cas_num*(1-Train)),1)*(-1);
for i=1:round(cas_num*(1-Train))
    for j =1:round(max(max(cascade_list))/60)
        if predictlabel2(i,j)==index(i+cas_num*Train,4)
            linear_time(i,1)=j*30;
            break;
        end
    end
end
for i=1:round(cas_num*(1-Train))
    if burst_time(i+cas_num*Train)~=-1
        if linear_time(i,1)==-1
            linear_time(i,1)=0 ;
        else
            linear_time(i,1) = max((burst_time(i+cas_num*Train)-linear_time(i,1)),0);
        end
    else
        linear_time(i,1)=-1;
    end
end

Prewhther_time = ones(round(cas_num*(1-Train)),1)*(-1);
for i=1:round(cas_num*(1-Train))
    for j =1:round(max(max(cascade_list))/60)
        if predictlabel3(i,j)==index(i+cas_num*Train,4)
            Prewhther_time(i,1)=j*30;
            break;
        end
    end
end
for i=1:round(cas_num*(1-Train))
    if burst_time(i+cas_num*Train)~=-1
        if Prewhther_time(i,1)==-1
            Prewhther_time(i,1)=0 ;
        else
            Prewhther_time(i,1) = max((burst_time(i+cas_num*Train)-Prewhther_time(i,1)),0);
        end
    else
        Prewhther_time(i,1)=-1;
    end
end

time1 = sum(svr_time(find(svr_time>0)))/length(find(svr_time>=0));
time2 = sum(linear_time(find(linear_time>0)))/length(find(linear_time>=0));
time3 = sum(Prewhther_time(find(Prewhther_time>0)))/length(find(Prewhther_time>=0));

%% ############################################ metric ############################################
for epoch = 1:6
    TP = 0;
    FP = 0;
    FN = 0;
    TN = 0;
    for i=1:cas_num*(1-Train)
        if predictlabel1(i,epoch)==1 & index(i+cas_num*Train,4)==1
            TP = TP + 1 ;
        else if predictlabel1(i,epoch)==1 & index(i+cas_num*Train,4)==0
                FP = FP + 1 ;
            else if predictlabel1(i,epoch)==0 & index(i+cas_num*Train,4)==1
                    FN = FN + 1 ;
                else
                    TN = TN +1 ;
                end
            end
        end
    end
    recall1(epoch) = TP/(TP+FN);
    precision1(epoch) = TP/(TP+FP);
    f11(epoch) = 2/(1/recall1(epoch)+1/precision1(epoch));
    [~,index_cover] = sort(index(1+cas_num*Train:end,3)-index(1+cas_num*Train:end,2));
    cover201(epoch) = sum(predictlabel1(index_cover(end-19:end),epoch))/20 ;
    auc1(epoch) = roc_curve(predictlabel1(:,epoch),index(:,4));
    TP = 0;
    FP = 0;
    FN = 0;
    TN = 0;
    for i=1:cas_num*(1-Train)
        if predictlabel2(i,epoch)==1 & index(i+cas_num*Train,4)==1
            TP = TP + 1 ;
        else if predictlabel2(i,epoch)==1 & index(i+cas_num*Train,4)==0
                FP = FP + 1 ;
            else if predictlabel2(i,epoch)==0 & index(i+cas_num*Train,4)==1
                    FN = FN + 1 ;
                else
                    TN = TN +1 ;
                end
            end
        end
    end
    recall2(epoch) = TP/(TP+FN);
    precision2 (epoch)= TP/(TP+FP);
    f12(epoch) = 2/(1/recall2(epoch)+1/precision2(epoch));
    cover202(epoch) = sum(predictlabel2(index_cover(end-19:end),epoch))/20 ;
    auc2(epoch) = roc_curve(predictlabel2(:,epoch),index(:,4));
    TP = 0;
    FP = 0;
    FN = 0;
    TN = 0;
    for i=1:cas_num*(1-Train)
        if predictlabel3(i,epoch)==1 & index(i+cas_num*Train,4)==1
            TP = TP + 1 ;
        else if predictlabel3(i,epoch)==1 & index(i+cas_num*Train,4)==0
                FP = FP + 1 ;
            else if predictlabel3(i,epoch)==0 & index(i+cas_num*Train,4)==1
                    FN = FN + 1 ;
                else
                    TN = TN +1 ;
                end
            end
        end
    end
    recall3(epoch) = TP/(TP+FN);
    precision3(epoch) = TP/(TP+FP);
    f13(epoch) = 2/(1/recall3(epoch)+1/precision3(epoch));
    [~,index_cover] = sort(index(1+cas_num*Train:end,3)-index(1+cas_num*Train:end,2));
    cover203(epoch) = sum(predictlabel3(index_cover(end-19:end),epoch))/20 ;
    auc3(epoch) = roc_curve(predictlabel3(:,epoch),index(:,4));
end


%% ################################# survival #############################################
interval = 60;
cumhaz = xlsread('cumhaz');
% plot(exp(-cumhaz))
b = [4.852*10^-2  8.195*10^-6 -5.629*10^-3] ;
ht = diff(cumhaz);

X1 = [] ;% retweet_num
X2 = [] ;% follower
X3 = [] ;% original timestamp
for i =1:cas_num
    for j=1:round(max(cascade_list(i,:))/interval)
        X1(i,j) = length(find(cascade_list(i,:)<=j*interval & cascade_list(i,:)>=0)) ;
        X2(i,j) = sum(follower_list(i,(find(cascade_list(i,:)<=j*interval & cascade_list(i,:)>=0 )))) ;
        X3(i,j) = index(i,1);
    end
end

% ######################## jifen ###################################
for i=1:cas_num
    for j = 1:round(max(max(cascade_list))/60)
        if X1(i,j)==0
            break
        end
    end
    if j~=1
        X1(i,j:round(max(max(cascade_list))/60))=X1(i,j-1)*ones(1,round(max(max(cascade_list))/60)-j+1);
    end
end
for i=1:cas_num
    for j = 1:round(max(max(cascade_list))/60)
        if X2(i,j)==0
            break
        end
    end
    if j~=1
        X2(i,j:round(max(max(cascade_list))/60))=X2(i,j-1)*ones(1,round(max(max(cascade_list))/60)-j+1);
    end
end
for i=1:cas_num
    for j = 1:round(max(max(cascade_list))/60)
        if X3(i,j)==0
            break
        end
    end
    if j~=1
        X3(i,j:round(max(max(cascade_list))/60))=X3(i,j-1)*ones(1,round(max(max(cascade_list))/60)-j+1);
    end
end

cumht(1:1000,1:4529)=0;
for i=1:1000
    for j=2:4529
        cumht(i,j) = cumht(i,j-1)+ht(j)*exp(b(1)*X1(i,floor(j/interval)+1)+b(2)*X2(i,floor(j/interval)+1)+b(3)*X3(i,floor(j/interval)+1));
    end
end
% ###########################################
cumht = exp(-cumht);
for i=1:4529
    tempV = cumht(find(index(1:800,4)==1 & cumht(1:800,2)~=0),i);
    tempN = cumht(find(index(1:800,4)==0 & cumht(1:800,2)~=0),i);
    [mu1,sigma1] = normfit(tempV);
    [mu2,sigma2] = normfit(tempN);
    S(i) = (mu1*max(sigma2,0.00001)+mu2*max(sigma1,0.00001))/(max(sigma1,0.00001)+max(sigma2,0.00001));
    Sansnet(i) = (length(find(index(1:800,4)==0))/800*sum(tempV)+length(find(index(1:800,4)==1))/800*sum(tempN))/(length(find(index(1:800,4)==0))/800*length(tempV)+length(find(index(1:800,4)==1))/800*length(tempN));
end
plot(S,'linewidth',2);
hold on ;
plot(cumht(find(index(801:1000,4)==0),:)','color','b')
plot(cumht(find(index(801:1000,4)==1),:)','color','r')
xlim([0 500]);
% 
% plot(Sansnet);


judgeS = zeros(6,200);
judgeSansnet = zeros(6,200);

% 
for i = find((index(801:1000,3)-index(801:1000,2))>=0.8*threshold)'+800
    for k =1:6
        for j =1:k*50
            if cumht(i,j)<S(j)
                judgeS(k,i-800) = 1;
            end
            if cumht(i,j)<Sansnet(j)
                judgeSansnet(k,i-800) = 1;
                break;
            end
        end
    end
end
judgeS = judgeS';
judgeSansnet = judgeSansnet';

for k=1:6
    [~,index_cover] = sort(index(801:1000,3)-index(801:1000,2));
    coverS(k) = sum(judgeS(k,index_cover(end-19:end)))/20 ;
    [~,index_cover] = sort(index(801:1000,3)-index(801:1000,2));
    coverSansnet(k) = sum(judgeSansnet(k,index_cover(end-19:end)))/20 ;
    recallS(k) = length(find(index(801:1000,4)==1 & judgeS(k,:)'==1))/length(find((index(801:1000,4)'==1)));
    recallSansnet(k) = length(find(index(801:1000,4)==1 & judgeSansnet(k,:)'==1))/length(find((index(801:1000,4)'==1)));
    precisionS(k) = length(find(index(801:1000,4)==1 & judgeS(k,:)'==1))/length(find((judgeS(k,:)'==1)));
    precisionSansnet(k) = length(find(index(801:1000,4)==1 & judgeSansnet(k,:)'==1))/length(find((judgeSansnet(k,:)'==1)));
    f1Sansnet(k) = 2/(1/recallSansnet(k)+1/precisionSansnet(k));
    f1S(k) = 2/(1/recallS(k)+1/precisionS(k));
    aucS(k) = roc_curve(judgeS(k,:)',index(:,4));
    aucSansnet(k) = roc_curve(judgeSansnet(k,:)',index(:,4));
end



%% ############################################# test ##############################################
result=[recall1; precision1; f11 ;cover201; auc1; recall2; precision2; f12 ;cover202; auc2; ...
    recall3 ;precision3 ;f13 ;cover203; auc3; recallSansnet; precisionSansnet; f1Sansnet ;coverSansnet; aucSansnet;...
    recallS; precisionS; f1S ;coverS; aucS];


result1 = [recall1(6)  precision1(6)  f11(6)  cover201(6)  auc1(6); recall2(6)  precision2(6) f12(6) cover202(6) auc2(6); ...
    recall3(6) precision3(6) f13(6) cover203(6) auc3(6); recallSansnet(6) precisionSansnet(6) f1Sansnet(6) coverSansnet(6) aucSansnet(6);...
    recallS(6) precisionS(6) f1S(6) coverS(6) aucS(6)];
