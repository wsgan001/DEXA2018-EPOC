X = zeros(cas_num,4);
time_window = 9;
hpara = [1.80559120130171;17.7423648387462;5.23202493942266;0.0193077598007735;3.69235695723155;0.872674014552439];
npara = [2.04814506271710;129.002787991689;2.83344931314487;0.0382115425063234;10.2529912426622;5.78913222646229] ;
for epoch=1:round(max(max(cascade_list))/60)
    fprintf('epoch %d start\n',epoch);
    observe_time = 60*epoch; %һСʱ
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
    po_c(:,1) = length(find(testX(1:cas_num*Train,4)==1))/(cas_num*Train) .* betapdf(testX(cas_num*Train+1:end,1),hpara(1),hpara(2)) .* gampdf(testX(cas_num*Train+1:end,2),hpara(3),hpara(4)) .* normpdf(testX(cas_num*Train+1:end,3),hpara(5),hpara(6));
    po_c(:,2) = length(find(testX(1:cas_num*Train,4)==0))/(cas_num*Train) .* betapdf(testX(cas_num*Train+1:end,1),npara(1),npara(2)) .* gampdf(testX(cas_num*Train+1:end,2),npara(3),npara(4)) .* normpdf(testX(cas_num*Train+1:end,3),npara(5),npara(6));
    pco = [po_c(:,1) ./ (po_c(:,1) + po_c(:,2)) , po_c(:,2) ./ (po_c(:,1) + po_c(:,2))];
    predictlabel3(1:round(cas_num*(1-Train)),epoch) = (pco(:,1)*0.9>=pco(:,2));
end


