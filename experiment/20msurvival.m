cumhaz = xlsread('cumhaz');
% plot(exp(-cumhaz))
b = [4.852*10^-2  8.195*10^-6 -5.629*10^-3] ;
ht = diff(cumhaz);

cumht(1:1000,1:4529)=0;
for i=1:1000
    for j=2:4529
        cumht(i,j) = cumht(i,j-1)+ht(j)*exp(b(1)*X1(i,floor(j/interval)+1)+b(2)*X2(i,floor(j/interval)+1)+b(3)*X3(i,floor(j/interval)+1));
    end
end

cumht = exp(-cumht);
% for i=1:4529
%     tempV = cumht(find(index(1:800,4)==1 & cumht(1:800,2)~=0),i);
%     tempN = cumht(find(index(1:800,4)==0 & cumht(1:800,2)~=0),i);
%     [mu1,sigma1] = normfit(tempV);
%     [mu2,sigma2] = normfit(tempN);
%     S(i) = (mu1*max(sigma2,0.00001)+mu2*max(sigma1,0.00001))/(max(sigma1,0.00001)+max(sigma2,0.00001));
%     Sansnet(i) = (length(find(index(1:800,4)==0))/800*sum(tempV)+length(find(index(1:8000,4)==1))/800*sum(tempN))/(length(find(index(1:800,4)==0))/800*length(tempV)+length(find(index(1:800,4)==1))/800*length(tempN));
% end
% plot(S,'linewidth',2);
% hold on ;
plot(cumht(find(index(10:50,4)==0),:)','color','b')
hold on;
plot(cumht(find(index(20:1000,4)==1),:)','color','r')
xlim([0 2500]);
% plot(Sansnet);


judgeS = zeros(1,2000);
judgeSansnet = zeros(1,2000);
for i = find((index(8001:10000,3)-index(8001:10000,2))>=0.8*threshold)'+8000
    for j =1:50
        if cumht(i,j)<S(j)
            judgeS(i-8000) = 1;
        end
        if cumht(i,j)<Sansnet(j)
            judgeSansnet(i-8000) = 1;
            break;
        end
    end
end

length(find(index(8001:10000,4)==1 & judgeS'==1))/length(find((index(8001:10000,4)'==1)))
length(find(index(8001:10000,4)==1 & judgeSansnet'==1))/length(find((index(8001:10000,4)'==1)))
length(find(index(8001:10000,4)==1 & judgeS'==1))/length(find((judgeS'==1)))
length(find(index(8001:10000,4)==1 & judgeSansnet'==1))/length(find((judgeSansnet'==1)))