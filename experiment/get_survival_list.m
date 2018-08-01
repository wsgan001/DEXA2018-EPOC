%% #############################################load data#######################################################
clear all 
clc
% data = textread('data.txt');
% index = textread('index.txt');
data = textread('data_small.txt');
index = textread('index_small.txt');
% index = index(1:1000,:);
% data = data(1:index(1000,3),:);
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
cascade_list = cascade_list/60 ; % ·ÖÖÓ
% follower_list = cumsum(follower_list')';

%% ############################################random sort############################################
shuffle = randperm(cas_num);
cascade_list = cascade_list(shuffle,:);
follower_list = follower_list(shuffle,:);
index = index(shuffle,:);
%% ############################################get threshold############################################
sum_rank = sort(index(1:cas_num*Train,3)-index(1:cas_num*Train,2)) ;
threshold = sum_rank(round(cas_num*Train*0.95));
index(:,4) = (index(:,3)-index(:,2))>threshold ;

%% ############################################ XTC ####################################################
interval = 60 ; %h
t = 1 ;
X = [] ;% cascade_num post_time total_follower
T = [] ;% tStart tStop
Censor = [] ;% 0/1
threshold2 = 50;
for i =1:cas_num*Train
    temp_list = cascade_list(i,cascade_list(i,:)>=0);
    temp_follower = follower_list(i,cascade_list(i,:)>=0);
    if temp_list(end)>interval
        if length(temp_list)>threshold2
            temp_list = temp_list(1:threshold2+1) ;
            temp_follower = follower_list(1:threshold2+1);
            if temp_list(end)>interval
                for j = 1:round(temp_list(end)/interval)
                    X(t,1) = length(temp_list( find((temp_list<j*interval))));
                    X(t,2) = index(i,1) ;
                    X(t,3) = sum(temp_follower( find((temp_list<j*interval))));
                    T(t,1) = (j-1)*interval ;
                    T(t,2) = j*interval ;
                    Censor(t) = 0 ;
                    t = t+1 ;
                end
                Censor(t-1) = 1 ;
                T(t-1,2) = temp_list(end);
            else
                X(t,1) = length(temp_list) ;
                X(t,2) = index(i,1) ;
                X(t,3) = sum(temp_follower);
                T(t,1) = 0 ;
                T(t,2) = temp_list(end) ;
                Censor(t) = 1 ;
                t = t+1 ;
            end
        else
            for j = 1:round(temp_list(end)/interval)
                X(t,1) = length(temp_list( find((temp_list<j*interval))));
                X(t,2) = index(i,1) ;
                X(t,3) = sum(temp_follower( find((temp_list<j*interval))));
                T(t,1) = (j-1)*interval ;
                T(t,2) = j*interval ;
                Censor(t) = 0 ;
                t = t+1 ;
            end
            T(t-1,2) = temp_list(end);
        end     
    else
        X(t,1) = length(temp_list) ;
        X(t,2) = index(i,1) ;
        X(t,3) = sum(temp_follower);
        T(t,1) = 0 ;
        T(t,2) = temp_list(end) ;
        if X(t,1)>threshold2
            Censor(t) = 1 ;
        else
            Censor(t) = 0 ;
        end
        t = t+1 ;
    end
end
Censor = Censor';
data = [X,T,Censor];