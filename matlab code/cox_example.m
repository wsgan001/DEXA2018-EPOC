load(fullfile(matlabroot,'examples','stats','lightbulb.mat'));
b = coxphfit(lightbulb(:,2),lightbulb(:,1), ...
'Censoring',lightbulb(:,3))
coxphopt = statset('coxphfit');
coxphopt.Display = 'final';
coxphopt.MaxIter = 50;
b = coxphfit(lightbulb(:,2),lightbulb(:,1),...
'Censoring',lightbulb(:,3),'Options',coxphopt)

rng('default') % for reproducibility
X = 4*rand(100,1);
A = 50*exp(-0.5*X);
B = 2;
y = wblrnd(A,B);
[b,logL,H,stats] = coxphfit(X,y);
[b logL]
stairs(H(:,1),exp(-H(:,2)),'LineWidth',2)