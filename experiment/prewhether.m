%AAAI2015
function  result = prewhether(trainX,trainY,testX,testY)

[Ntrain,~] = size(trainX);
[Ntest,~] = size(testX);
pc = [sum(trainY)/Ntrain,1-sum(trainY)/Ntrain];
po_c = zeros(Ntest,2);
pco = zeros(Ntest,2);

options = optimoptions('fmincon','Display','off');
[hpara,fval] = fmincon(@(para) likelihood_b1_1(trainX(find(trainY==1)),para),rand(6,1),[],[],[],[],[0,0,0,0,-200,0],[200,200,200,200,200,200], [] ,options);
[npara,fval] = fmincon(@(para) likelihood_b1_1(trainX(find(trainY==0)),para),rand(6,1),[],[],[],[],[0,0,0,0,-200,0],[200,200,200,200,200,200], [],options);

po_c(:,1) = pc(1) .* betapdf(testX(:,1),hpara(1),hpara(2)) .* gampdf(testX(:,2),hpara(3),hpara(4)) .* normpdf(testX(:,3),hpara(5),hpara(6));
po_c(:,2) = pc(2) .* betapdf(testX(:,1),npara(1),npara(2)) .* gampdf(testX(:,2),npara(3),npara(4)) .* normpdf(testX(:,3),npara(5),npara(6));
pco = [po_c(:,1) ./ (po_c(:,1) + po_c(:,2)) , po_c(:,2) ./ (po_c(:,1) + po_c(:,2))];
predictlabel = (pco(:,1)>=pco(:,2));

%Precision/Recall/F1 score
precision = sum((predictlabel==1)&(testY==1))/sum(predictlabel==1);
recall = sum((predictlabel==1)&(testY==1))/sum(testY==1);
F1 = 2/(1/precision+1/recall);
% auc = roc_curve(pco(:,1),testY(:,2));

%hot topic recall
w = 20;
[~,location] = sort(testY(:,1),'descend');
topk_number = location(1:w);
recall1 = sum((predictlabel(topk_number)==1)&(testY(topk_number,2)==1))/sum(testY(topk_number,2)==1);
%auc1 = roc_curve(pco(topk_number,1),testY(topk_number,2));

result = [precision, recall, accuracy, F1, auc, precision1, recall1, accuracy1, F11, mse, mapd, mpd, md, nse, pearson1, spearman1, kendall1, pearson2, spearman2, kendall2];

end