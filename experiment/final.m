%% #############################################load data#######################################################
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
%% ######################################################list cascade############################################
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
% shuffle = randperm(cas_num);
% cascade_list = cascade_list(shuffle,:);
% follower_list = follower_list(shuffle,:);
% index = index(shuffle,:);
%% ############################################get threshold############################################
sum_rank = sort(index(1:cas_num*Train,3)-index(1:cas_num*Train,2)) ;
threshold = sum_rank(round(cas_num*Train*0.95));
index(:,4) = (index(:,3)-index(:,2))>threshold ;
% t = 1 ;
% for i =1:cas_num
%     temp = cascade_list(i,threshold);
%     if temp>0
%         burst_time(t) = temp;
%         original_time(t) = index(i,1);
%         cascade_list2(t,:) = cascade_list(i,:);
%         follower_list2(t,:) = follower_list(i,:);
%         t = t+1;
%     else
% end
for i =1:cas_num
    temp = cascade_list(i,threshold);
    if temp>0
        burst_time(i) = temp;
    else
        burst_time(i) = -1;
    end
end

%% ############################################get feature: total_sum, total_folower, orignial time, Ylabel############################################
X = zeros(cas_num,4);
time_window = 9;
% hpara = [1.80559120130171;17.7423648387462;5.23202493942266;0.0193077598007735;3.69235695723155;0.872674014552439];
% npara = [2.04814506271710;129.002787991689;2.83344931314487;0.0382115425063234;10.2529912426622;5.78913222646229] ;
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
    options = optimoptions('fmincon','Display','off');
    if epoch==1
        [hpara,~] = fmincon(@(para) likelihood_b1_1(testX(find(testX(1:cas_num*Train,4)==1),1:3),para),rand(6,1),[],[],[],[],[0,0,0,0,-200,0],[200,200,200,200,200,200], [] ,options);
        [npara,~] = fmincon(@(para) likelihood_b1_1(testX(find(testX(1:cas_num*Train,4)==0),1:3),para),rand(6,1),[],[],[],[],[0,0,0,0,-200,0],[200,200,200,200,200,200], [],options);
        fprintf('PreWhether finished!\n');
    end
    po_c(:,1) = length(find(testX(1:cas_num*Train,4)==1))/(cas_num*Train) .* betapdf(testX(cas_num*Train+1:end,1),hpara(1),hpara(2)) .* gampdf(testX(cas_num*Train+1:end,2),hpara(3),hpara(4)) .* normpdf(testX(cas_num*Train+1:end,3),hpara(5),hpara(6));
    po_c(:,2) = length(find(testX(1:cas_num*Train,4)==0))/(cas_num*Train) .* betapdf(testX(cas_num*Train+1:end,1),npara(1),npara(2)) .* gampdf(testX(cas_num*Train+1:end,2),npara(3),npara(4)) .* normpdf(testX(cas_num*Train+1:end,3),npara(5),npara(6));
    pco = [po_c(:,1) ./ (po_c(:,1) + po_c(:,2)) , po_c(:,2) ./ (po_c(:,1) + po_c(:,2))];
    predictlabel3(1:round(cas_num*(1-Train)),epoch) = (pco(:,1)*0.9>=pco(:,2));
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
            linear_time(i,1)=j*60;
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
            Prewhther_time(i,1)=j*60;
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

epoch = 1;
%% ############################################ metric############################################
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
recall1 = TP/(TP+FN);
precision1 = TP/(TP+FP);
f11 = 2/(1/recall1+1/precision1);
[~,index_cover] = sort(index(1+cas_num*Train:end,3)-index(1+cas_num*Train:end,2));
cover201 = sum(predictlabel1(index_cover(end-19:end),epoch))/20 ;
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
recall2 = TP/(TP+FN);
precision2 = TP/(TP+FP);
f12 = 2/(1/recall2+1/precision2);
cover202 = sum(predictlabel2(index_cover(end-19:end),epoch))/20 ;
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
recall3 = TP/(TP+FN);
precision3 = TP/(TP+FP);
f13 = 2/(1/recall3+1/precision3);
[~,index_cover] = sort(index(1+cas_num*Train:end,3)-index(1+cas_num*Train:end,2));
cover203 = sum(predictlabel3(index_cover(end-19:end),epoch))/20 ;
% fprintf('Linear classification finish.\n')
% %% get feature: total_sum, total_folower, orignial time, burst_time
% X2 = zeros(t-1,4);
% observe_time = 24*60; %一天
% for i =1:t-1
%     X2(i,1) = length(find(cascade_list2(i,:)<observe_time));
%     X2(i,2) = sum(follower_list2(find(cascade_list2(i,:)<observe_time)));
%     X2(i,3) = original_time(i);
%     X2(i,4) = burst_time(i);
% end
% 
% %% SVR regression
% baseline3 = fitrsvm(X2(1:round((t-1)*Train),1:3),X2(1:round((t-1)*Train),4),'KernelFunction','gaussian','Standardize',true);
% predict_time1 = predict(baseline3,X2(round((t-1)*Train)+1:end,1:3));
% predict1 = max((burst_time(round((t-1)*Train)+1:end)'-predict_time1),0);
% time1 = sum(predict1)/length(predict1);
% fprintf('SVR regression finish.\n')
% %% Linear regression
% baseline4 = fitlm(X2(1:round((t-1)*Train),1:3),X2(1:round((t-1)*Train),4));
% predict_time2 = predict(baseline4,X2(round((t-1)*Train)+1:end,1:3));
% predict2 = max((burst_time(round((t-1)*Train)+1:end)'-predict_time2),0);
% time2 = sum(predict2)/length(predict2);
% fprintf('Linear regression finish.\n')
% 



%% ############################################get XTC list############################################
% interval = 20 ; %h
% t = 1 ;
% X = [] ;% cascade_num post_time total_follower
% T = [] ;% tStart tStop
% Censor = [] ;% 0/1
% threshold2 = 50;
% for i =1:cas_num*Train
%     temp_list = cascade_list(i,cascade_list(i,:)>=0);
%     temp_follower = follower_list(i,cascade_list(i,:)>=0);
%     if temp_list(end)>interval
%         if length(temp_list)>threshold2
%             temp_list = temp_list(1:threshold2+1) ;
%             temp_follower = follower_list(1:threshold2+1);
%             if temp_list(end)>interval
%                 for j = 1:round(temp_list(end)/interval)
%                     X(t,1) = length(temp_list( find((temp_list<j*interval))));
%                     X(t,2) = index(i,1) ;
%                     X(t,3) = sum(temp_follower( find((temp_list<j*interval))));
%                     T(t,1) = (j-1)*interval ;
%                     T(t,2) = j*interval ;
%                     Censor(t) = 0 ;
%                     t = t+1 ;
%                 end
%                 Censor(t-1) = 1 ;
%                 T(t-1,2) = temp_list(end);
%             else
%                 X(t,1) = length(temp_list) ;
%                 X(t,2) = index(i,1) ;
%                 X(t,3) = sum(temp_follower);
%                 T(t,1) = 0 ;
%                 T(t,2) = temp_list(end) ;
%                 Censor(t) = 1 ;
%                 t = t+1 ;
%             end
%         else
%             for j = 1:round(temp_list(end)/interval)
%                 X(t,1) = length(temp_list( find((temp_list<j*interval))));
%                 X(t,2) = index(i,1) ;
%                 X(t,3) = sum(temp_follower( find((temp_list<j*interval))));
%                 T(t,1) = (j-1)*interval ;
%                 T(t,2) = j*interval ;
%                 Censor(t) = 0 ;
%                 t = t+1 ;
%             end
%             T(t-1,2) = temp_list(end);
%         end     
%     else
%         X(t,1) = length(temp_list) ;
%         X(t,2) = index(i,1) ;
%         X(t,3) = sum(temp_follower);
%         T(t,1) = 0 ;
%         T(t,2) = temp_list(end) ;
%         if X(t,1)>threshold2
%             Censor(t) = 1 ;
%         else
%             Censor(t) = 0 ;
%         end
%         t = t+1 ;
%     end
% end
% Censor = Censor';
% data = [X,T,Censor];
% 
% %% ############################################# survival boundary ###################################################
X1 = [] ;% cascade_num
X2 = [] ;% follower
X3 = [] ;% original timestamp
interval = 60;
for i =1:cas_num
    for j=1:round(max(max(cascade_list))/interval)
        X1(i,j) = length(find(cascade_list(i,:)<=j*interval & cascade_list(i,:)>=0)) ;
        X2(i,j) = sum(follower_list(i,(find(cascade_list(i,:)<=j*interval & cascade_list(i,:)>=0 )))) ;
        X3(i,j) = index(i,1);
    end
    fprintf('cascade %d finish\n',i);
end
b = [4.852*10^-2  8.195*10^-6 -5.629*10^-3] ;
% % Separation= []; % survival boundary
% for t = 1:round(max(max(cascade_list))/interval)
%     tempV1 = X1(find(index(1:cas_num*Train,4)==1),t);
%     tempN1 = X1(find(index(1:cas_num*Train,4)==0),t);
%     tempV2 = X2(find(index(1:cas_num*Train,4)==1),t);
%     tempN2 = X2(find(index(1:cas_num*Train,4)==0),t);
%     tempV3 = X3(find(index(1:cas_num*Train,4)==1),t);
%     tempN3 = X3(find(index(1:cas_num*Train,4)==0),t);
%     tempV4 = [tempV1 tempV2 tempV3];
%     tempN4 = [tempN1 tempN2 tempN3];
%     tempV(:,t) = exp(-exp(tempV4*(b')));
%     tempN(:,t) = exp(-exp(tempN4*(b')));
% end
% A = zeros(round(max(max(cascade_list))/interval)-1,round(max(max(cascade_list))/interval));
% for i = 1:round(max(max(cascade_list))/interval)-1
%    A(i,i) = -1 ;
%    A(i,i+1) = 1;
% end
% options = optimoptions('fmincon','Display','off');
% [Separation,~] = fmincon(@(para) separation(tempV,tempN,para,round(max(max(cascade_list))/interval)),...
%     linspace(0.03,0.001,round(max(max(cascade_list))/interval)),A,zeros(round(max(max(cascade_list))/interval)-1,1),...
%     [],[],zeros(round(max(max(cascade_list))/interval),1),...
%     ones(round(max(max(cascade_list))/interval),1),[],options);


S= []; % survival boundary
for t = 1:round(max(max(cascade_list))/interval)
    tempV1 = X1(find(index(1:cas_num*Train,4)==1),t);
    tempN1 = X1(find(index(1:cas_num*Train,4)==0),t);
    tempV2 = X2(find(index(1:cas_num*Train,4)==1),t);
    tempN2 = X2(find(index(1:cas_num*Train,4)==0),t);
    tempV3 = X3(find(index(1:cas_num*Train,4)==1),t);
    tempN3 = X3(find(index(1:cas_num*Train,4)==0),t);
    tempV = [tempV1 tempV2 tempV3];
    tempN = [tempN1 tempN2 tempN3];
    tempV = exp(-exp(tempV*(b')));
    tempN = exp(-exp(tempN*(b')));
    [muV,sigmaV] = normfit(tempV);
    [muN,sigmaN] = normfit(tempN);
    S(t) = (muV*sigmaN+muN*sigmaV)/(sigmaV+sigmaN) ;
    Separation(t) = (sum(tempV)/length(tempV)*length(tempN)+sum(tempN)/length(tempN)*length(tempV))/(length(tempV)+length(tempN));
end
%% ############################################# test ##############################################

% ###########################Sansnet############################
temp=[] ;
model_sansnet = ones(round(cas_num,1),1)*(-1);
for i = cas_num*Train+1:cas_num
    for t = 1:round(max(max(cascade_list))/interval)
        temp(t,1) = X1(i,t);
        temp(t,2) = X2(i,t);
        temp(t,3) = X3(i,t);
    end
    fitting = exp(-exp(temp*b')) ;
    for t = 1:round(max(max(cascade_list))/interval)
        if fitting(t)<Separation(t)
            model_sansnet(i) = t;
            break
        end
    end
end

model_sansnet_time = ones(round(cas_num,1),1)*(-1) ;
for i = cas_num*Train+1:cas_num
    if index(i,4)==1 & model_sansnet(i)>0
        if burst_time(i) - model_sansnet(i)*interval>0
            model_sansnet_time(i) = burst_time(i) - model_sansnet(i)*interval;
        else
            model_sansnet_time(i) = 0;
        end
    else if index(i,4)==1 & model_sansnet(i)==-1
            model_sansnet_time(i) = 0 ;
        else 
            model_sansnet_time(i) = -1 ;
        end
    end
end

model_sansnet_time = max(model_sansnet_time(find(model_sansnet_time>=0)),0);
time4 = sum(model_sansnet_time)/length(model_sansnet_time);
predictlabel4 = (model_sansnet(cas_num*Train+1:cas_num)==3 | model_sansnet(cas_num*Train+1:cas_num)==1 | model_sansnet(cas_num*Train+1:cas_num)==2);
TP = 0;
FP = 0;
FN = 0;
TN = 0;
for i=1:cas_num*(1-Train)
    if predictlabel4(i)==1 & index(i+cas_num*Train,4)==1
        TP = TP + 1 ;
    else if predictlabel4(i)==1 & index(i+cas_num*Train,4)==0
            FP = FP + 1 ;
        else if predictlabel4(i)==0 & index(i+cas_num*Train,4)==1
                FN = FN + 1 ;
            else
                TN = TN +1 ;
            end
        end
    end
end
recall4 = TP/(TP+FN);
precision4 = TP/(TP+FP);
f14 = 2/(1/recall4+1/precision4);
[~,index_cover] = sort(index(1+cas_num*Train:end,3)-index(1+cas_num*Train:end,2));
cover204 = sum(predictlabel4(index_cover(end-19:end),epoch))/20 ;



% ########################## EPOC ####################################
threshold2 = 1;
% boundary
temp=[] ;
model = ones(round(cas_num,1),1)*(-1);
for i = cas_num*Train+1:cas_num
    for t = 1:round(max(max(cascade_list))/interval)
        temp(t,1) = X1(i,t);
        temp(t,2) = X2(i,t);
        temp(t,3) = X3(i,t);
    end
    fitting = exp(-exp(temp*b')) ;
    for t = 1:round(max(max(cascade_list))/interval)
        if fitting(t)<S(t)
            model(i) = t;
            break
        end
    end
end

model_time = ones(round(cas_num,1),1)*(-1) ;
for i = cas_num*Train+1:cas_num
    if index(i,4)==1 & model(i)>0
        if burst_time(i) - model(i)*interval>0
            model_time(i) = burst_time(i) - model(i)*interval;
        else
            model_time(i) = 0;
        end
    else if index(i,4)==1 & model(i)==-1
            model_time(i) = 0 ;
        else 
            model_time(i) = -1 ;
        end
    end
end

% hazard ceiling
temp=[] ;
model2 = ones(round(cas_num,1),1)*(-1);
for i = cas_num*Train+1:cas_num
    for t = 1:round(max(max(cascade_list))/interval)
        temp(t,1) = X1(i,t);
        temp(t,2) = X2(i,t);
        temp(t,3) = X3(i,t);
    end
    fitting = exp(-exp(temp*b')) ;
    hazard = -[0; diff(fitting)]./fitting ;
    for t = 1:round(max(max(cascade_list))/interval)
        if hazard(t)>threshold2
            model2(i) = t;
            break
        end
    end
end

model_time2 = ones(round(cas_num,1),1)*(-1) ;
for i = cas_num*Train+1:cas_num
    if index(i,4)==1 & model2(i)>0
        if burst_time(i) - model2(i)*interval>0
            model_time2(i) = burst_time(i) - model2(i)*interval;
        else
            model_time2(i) = 0;
        end
    else if index(i,4)==1 & model2(i)==-1
            model_time2(i) = 0 ;
        else 
            model_time2(i) = -1 ;
        end
    end
end


model_time = max(model_time(find(model_time>=0)),0);
model_time2 = max(model_time2(find(model_time2>=0)),0);
model_time = max(model_time,model_time2) ;

time5 = sum(model_time)/length(model_time);
predictlabel5 = (model(cas_num*Train+1:cas_num)==3 | model(cas_num*Train+1:cas_num)==1 | model(cas_num*Train+1:cas_num)==2);
TP = 0;
FP = 0;
FN = 0;
TN = 0;
for i=1:cas_num*(1-Train)
    if predictlabel4(i)==1 & index(i+cas_num*Train,4)==1
        TP = TP + 1 ;
    else if predictlabel4(i)==1 & index(i+cas_num*Train,4)==0
            FP = FP + 1 ;
        else if predictlabel4(i)==0 & index(i+cas_num*Train,4)==1
                FN = FN + 1 ;
            else
                TN = TN +1 ;
            end
        end
    end
end
recall5 = TP/(TP+FN);
precision5 = TP/(TP+FP);
f15 = 2/(1/recall5+1/precision5);
[~,index_cover] = sort(index(1+cas_num*Train:end,3)-index(1+cas_num*Train:end,2));
cover205 = sum(predictlabel5(index_cover(end-19:end),epoch))/20 ;
result=[recall1 precision1 f11 cover201 time1; recall2 precision2 f12 cover202 time2;...
    recall3 precision3 f13 cover203 time3; recall4 precision4 f14 cover204 time4;...
    recall5 precision5 f15 cover205 time5];

