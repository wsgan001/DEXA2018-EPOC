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
%% ############################################random sort############################################
shuffle = randperm(cas_num);
cascade_list = cascade_list(shuffle,:);
follower_list = follower_list(shuffle,:);
index = index(shuffle,:);
%% ############################################get threshold############################################
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

%% ############################################get feature: total_sum, total_folower, orignial time, Ylabel############################################
X = zeros(cas_num,4);
time_window = 9;
hpara = [0.561350005235447;3.72555275765086;1.01799621925657;0.218723556308548;1.47073250108350;1.47569893481172];
npara = [1.30622245904933;36.0318215995954;0.856594179736867;0.423740497546436;10.1176126311958;6.99545513382115] ;
for epoch=1:20
    fprintf('epoch %d start\n',epoch);
    observe_time = 30*epoch; %一小时
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
    predict_count1 = predict(baseline1,X(:,1:3));
    predictlabel1(:,epoch) = (predict_count1>=threshold);
    % Linear
    baseline2 = fitlm(X(1:cas_num*Train,1:3),X(1:cas_num*Train,4));
    predict_count2 = predict(baseline2,X(:,1:3));
    predictlabel2(:,epoch) = (predict_count2>=threshold);
    % Prewhether
    options = optimoptions('fmincon','Display','off');
    if 1==2
        [hpara,~] = fmincon(@(para) likelihood_b1_1(testX(find(testX(1:cas_num*Train,4)==1),1:3),para),rand(6,1),[],[],[],[],[0,0,0,0,-200,0],[200,200,200,200,200,200], [] ,options);
        [npara,~] = fmincon(@(para) likelihood_b1_1(testX(find(testX(1:cas_num*Train,4)==0),1:3),para),rand(6,1),[],[],[],[],[0,0,0,0,-200,0],[200,200,200,200,200,200], [],options);
        fprintf('PreWhether finished!\n');
    end
    po_c(:,1) = length(find(testX(1:cas_num*Train,4)==1))/(cas_num*Train) .* betapdf(testX(:,1),hpara(1),hpara(2)) .* gampdf(testX(:,2),hpara(3),hpara(4)) .* normpdf(testX(:,3),hpara(5),hpara(6));
    po_c(:,2) = length(find(testX(1:cas_num*Train,4)==0))/(cas_num*Train) .* betapdf(testX(:,1),npara(1),npara(2)) .* gampdf(testX(:,2),npara(3),npara(4)) .* normpdf(testX(:,3),npara(5),npara(6));
    pco = [po_c(:,1) ./ (po_c(:,1) + po_c(:,2)) , po_c(:,2) ./ (po_c(:,1) + po_c(:,2))];
    predictlabel3(:,epoch) = (pco(:,1)>=pco(:,2));
end


svr_time = ones(cas_num,1)*(-1);
for i=1:cas_num
    for j =1:12
        if predictlabel1(i,j)==index(i,4)
            svr_time(i,1)=j*60;
            break;
        end
    end
end
for i=1:cas_num
    if burst_time(i)~=-1
        if svr_time(i,1)==-1
            svr_time(i,1)=0 ;
        else
            svr_time(i,1) = max((burst_time(i)-svr_time(i,1))/burst_time(i),0);
        end
    else
        svr_time(i,1)=-1;
    end
end

linear_time = ones(cas_num,1)*(-1);
for i=1:cas_num
    for j =1:12
        if predictlabel2(i,j)==index(i,4)
            linear_time(i,1)=j*60;
            break;
        end
    end
end
for i=1:cas_num
    if burst_time(i)~=-1
        if linear_time(i,1)==-1
            linear_time(i,1)=0 ;
        else
            linear_time(i,1) = max((burst_time(i)-linear_time(i,1))/burst_time(i),0);
        end
    else
        linear_time(i,1)=-1;
    end
end

Prewhther_time = ones(cas_num,1)*(-1);
for i=1:cas_num
    for j =1:12
        if predictlabel3(i,j)==index(i,4)
            Prewhther_time(i,1)=j*60;
            break;
        end
    end
end
for i=1:cas_num
    if burst_time(i)~=-1
        if Prewhther_time(i,1)==-1
            Prewhther_time(i,1)=0 ;
        else
            Prewhther_time(i,1) = max((burst_time(i)-Prewhther_time(i,1))/burst_time(i),0);
        end
    else
        Prewhther_time(i,1)=-1;
    end
end

time1 = sum(svr_time(find(svr_time>0)))/length(find(svr_time>=0));
time2 = sum(linear_time(find(linear_time>0)))/length(find(linear_time>=0));
time3 = sum(Prewhther_time(find(Prewhther_time>0)))/length(find(Prewhther_time>=0));

%% ############################################ metric############################################
cost1 = 5;
cost2 = 1;
for epoch = 1:20
TP = 0;
FP = 0;
FN = 0;
TN = 0;
for i=1:cas_num
    if predictlabel1(i,epoch)==1 & index(i,4)==1
        TP = TP + 1 ;
    else if predictlabel1(i,epoch)==1 & index(i,4)==0
            FP = FP + 1 ;
        else if predictlabel1(i,epoch)==0 & index(i,4)==1
                FN = FN + 1 ;
            else
                TN = TN +1 ;
            end
        end
    end
end
% recall1 = TP/(TP+FN);
% precision1 = TP/(TP+FP);
% f11 = 2/(1/recall1+1/precision1);
% [~,index_cover] = sort(index(1:end,3)-index(1:end,2));
% cover201 = sum(predictlabel1(index_cover(end-19:end),epoch))/20 ;
FNR = FN/(FN+TP);
FPR = FP/(FP+TN);
p = (TP+FN)/1000 ;
cost11 = (FNR*p*cost1+FPR*(1-p)*cost2)/(p*cost1+(1-p)*cost2);
TP = 0;
FP = 0;
FN = 0;
TN = 0;
for i=1:cas_num
    if predictlabel2(i,epoch)==1 & index(i,4)==1
        TP = TP + 1 ;
    else if predictlabel2(i,epoch)==1 & index(i,4)==0
            FP = FP + 1 ;
        else if predictlabel2(i,epoch)==0 & index(i,4)==1
                FN = FN + 1 ;
            else
                TN = TN +1 ;
            end
        end
    end
end
% recall2 = TP/(TP+FN);
% precision2 = TP/(TP+FP);
% f12 = 2/(1/recall2+1/precision2);
% cover202 = sum(predictlabel2(index_cover(end-19:end),epoch))/20 ;
FNR = FN/(FN+TP);
FPR = FP/(FP+TN);
p = (TP+FN)/1000 ;
cost12 = (FNR*p*cost1+FPR*(1-p)*cost2)/(p*cost1+(1-p)*cost2);
TP = 0;
FP = 0;
FN = 0;
TN = 0;
for i=1:cas_num
    if predictlabel3(i,epoch)==1 & index(i,4)==1
        TP = TP + 1 ;
    else if predictlabel3(i,epoch)==1 & index(i,4)==0
            FP = FP + 1 ;
        else if predictlabel3(i,epoch)==0 & index(i,4)==1
                FN = FN + 1 ;
            else
                TN = TN +1 ;
            end
        end
    end
end
% recall3 = TP/(TP+FN);
% precision3 = TP/(TP+FP);
% f13 = 2/(1/recall3+1/precision3);
% [~,index_cover] = sort(index(1:end,3)-index(1:end,2));
% cover203 = sum(predictlabel3(index_cover(end-19:end),epoch))/20 ;
FNR = FN/(FN+TP);
FPR = FP/(FP+TN);
p = (TP+FN)/1000 ;
cost13 = (FNR*p*cost1+FPR*(1-p)*cost2)/(p*cost1+(1-p)*cost2);


% result(epoch,:)=[recall1 precision1 f11 cover201 time1 recall2 precision2 f12 cover202 time2 ...
%     recall3 precision3 f13 cover203 time3];
result(epoch,:)=[cost11 cost12 cost13];
end