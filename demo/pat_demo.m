% Demonstration of photometry regression analysis

% event onsets
data(1).events(1).name = 'cue';
data(1).events(1).ons = [1 10];
data(1).events(1).pmod(1).name = 'value';
data(1).events(1).pmod(1).param = [0.3 -0.9];

data(1).events(2).name = 'outcome';
data(1).events(2).ons = [3 13];
data(1).events(2).pmod(1).name = 'RPE';
data(1).events(2).pmod(1).param = [0.1 0.5];

% design matrix
[X, name] = pat_design(data);

% construct hypothetical neural data
b = [0 25 100 50 80]';
sd = 0.001;
y = normrnd(X*b,sd);

% run regression to recover coefficients
results = pat_regress(y,X);

% plot results
figure;
subplot(1,2,1);
plot(b,results.b,'+k','LineWidth',4,'MarkerSize',10);
xlabel('True coefficients','FontSize',25);
ylabel('Estimated coefficients','FontSize',25);
set(gca,'FontSize',20,'XLim',[-10 110],'YLim',[-10 110]);
axis square
subplot(1,2,2);
plot(y,'-r','LineWidth',4); hold on; plot(results.yhat,'-k','LineWidth',4)
xlabel('Time','FontSize',25);
ylabel('Calcium signal','FontSize',25);
legend({'Data' 'Model'}','FontSize',25);
set(gca,'FontSize',20);
set(gcf,'Position',[200 200 1200 500]);