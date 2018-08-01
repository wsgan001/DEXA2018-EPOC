clear all
clc
% https://www.mathworks.com/help/stats/cox-proportional-hazards-model-with-time-dependent-covariates.html?searchHighlight=cox&s_tid=doc_srchtitle
% https://www.mathworks.com/help/stats/readmission-times.html
load(fullfile(matlabroot,'examples','stats','simplesurvivaldata.mat'));
load(fullfile(matlabroot,'examples','stats','survivaldatacp.mat'));


mTime = [0 50 100]; % Measurement time
threeLabs = [labS.Lab_0 labS.Lab_50 labS.Lab_100];
nLabMeasure = sum(sum(~isnan(threeLabs))); % Number of lab measurements
data = zeros(nLabMeasure,6); % One row for each observation
oID = 0; % Observation ID
for i = 1 : size(labS,1)
    idx = find(mTime <= labS.Time(i));
    for j = 1 : length(idx)-1
        oID = oID + 1;
        data(oID,:) = [labS.ID(i) mTime(j:j+1) 1 labS.Sex(i) threeLabs(i,j)];
    end
    oID = oID + 1;
    data(oID,:) = [labS.ID(i) mTime(length(idx)) labS.Time(i) ...
            labS.Censoring(i) labS.Sex(i) threeLabs(i,length(idx))];
end
labCP = table(data(:,1),data(:,2),data(:,3),data(:,4),data(:,5),data(:,6), ...
    'VariableNames', {'ID','tStart','tStop','Censoring','Sex','Lab'});
idxInvalid = labCP.ID(find(labCP.tStart == labCP.tStop)) ;
idxAdjust = find(labCP.ID==idxInvalid);
labCP.tStop(idxAdjust(1)) = labCP.tStop(idxAdjust(1))-0.5;
labCP.tStart(idxAdjust(2)) = labCP.tStart(idxAdjust(2))-0.5;
X = [labCP.Sex labCP.Lab];
T = [labCP.tStart labCP.tStop];
[b,logL,H,stats] = coxphfit(X,T,'Censoring',labCP.Censoring,'Baseline',0)

stairs(H(:,1),exp(-H(:,2)),'LineWidth',2)
