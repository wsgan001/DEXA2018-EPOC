function y=separation(tempV, tempN, para,S)
    [V,~] = size(tempV);
    [N,~] = size(tempN);
    yhat = 0 ;
    for j = 1:S
        for i=1:V
            yhat = yhat+N/(N+V)*abs(tempV(i,j)-para(j))/(tempV(i,j)-para(j)) ;
        end
        for i=1:N
            yhat = yhat+V/(N+V)*abs(tempN(i,j)-para(j))/(tempN(i,j)-para(j))*(-1) ;
        end
    end
    y = yhat;
end