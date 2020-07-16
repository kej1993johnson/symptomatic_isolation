% Use this script to digitize image from paper where 2sd around measurement
% of some data ( in this case diameter of a CT scan) is measured

figflu = imread('Ip_Fig2.jpg');
figSARS = imread('Peires_Fig4.jpg');
% digitize 2 returns the x values (first column) and y values (second
% column) of the points that you have acquired
%% SARS dynamics
figdigSARS = digitize2('Peires_Fig4.jpg');
save('../out/SARS.mat', 'figdigSARS');
%% ARI flu dynamics figure
figdigARI = digitize2('Ip_Fig2.jpg');
save('../out/fluARI.mat', 'figdig1');
%% Pauci symptomatic flu dynamics
figdigpauci = digitize2('Ip_Fig2.jpg');
save('../out/flupauci.mat','figdigpauci');
%% He et al SARS
figdigHeSARS = digitize2('Hefig1A.jpg');
save('../out/SARSinf_t.mat', 'figdigHeSARS')
%% He et al flu
figdigHeflu = digitize2('Hefig1A.jpg');
save('../out/fluinf_t.mat', 'figdigHeflu')
%% Plot the digitized data
figure;
subplot(1,2,1)
plot(figdigHeSARS(:,1), figdigHeSARS(:,2),'*-', 'LineWidth', 2)
hold on
xlabel('days since symptom onset')
ylabel('viral load (RT-PCR copies/mL)')
legend(' SARS He et al', 'flu ARI Ip et al')
legend boxoff
title('SARS infectiousness (not normalized)')
set(gca,'FontSize',16,'LineWidth',1.5)
ylim([ 0 0.6])
subplot(1,2,2)
plot(figdigHeflu(:,1), figdigHeflu(:,2),'*-', 'LineWidth', 2)
hold on
xlabel('days since symptom onset')
ylabel('viral load (RT-PCR copies/mL)')
legend(' Flu He et al', 'flu ARI Ip et al')
legend boxoff
title('SARS infectiousness (not normalized)')
set(gca,'FontSize',16,'LineWidth',1.5)
ylim([ 0 0.6])


%% Fit flu to normal and gamma distributions
disease = 'flu';
distrib = 'gamma';
ydata = figdigHeflu(:,2);
ydata = ydata./sum(ydata);
ydata = ydata- min(ydata);
torig = figdigHeflu(:,1);
mint = min(torig);
t = torig-mint;
p0 = [3.5, 1];
yguess = gampdf(t, p0(1), p0(2)); 

% Fit flu to gamma 
err = @(p) objfun(p,t, ydata,disease, distrib);
[pbest, resnorm, residuals] = lsqnonlin(err, p0);
tmod = t(1):0.1:t(end);
ymod = gampdf(tmod, pbest(1), pbest(2));
sumsqerr = sum((residuals.^2));
figure;
subplot(1,2,1)
plot(tmod+mint, ymod, '-', 'LineWidth', 2)
hold on
plot(torig, ydata, 'o', 'LineWidth', 2)
legend('gamma model fit', 'data')
legend boxoff
set(gca,'FontSize',16,'LineWidth',1.5)
xlabel ('days since symptom onset')
ylabel('Density')
title(['Flu gamma fit sum-squared error =', num2str(sumsqerr)])
save('../out/flu_gammaparams.mat', 'pbest')
%
disease = 'flu';
distrib = 'normal';
ydata = figdigHeflu(:,2);
ydata = ydata./sum(ydata);
ydata = ydata- min(ydata);
torig = figdigHeflu(:,1);
mint = min(torig);
t = torig-mint;
p0 = [3.5, 1];
yguess = normpdf(t, p0(1), p0(2)); 

% Fit flu to normal
err = @(p) objfun(p,t, ydata,disease, distrib);
[pbestnorm, resnorm, residuals] = lsqnonlin(err, p0);
tmod = t(1):0.1:t(end);
ymod = normpdf(tmod, pbest(1), pbest(2));
sumsqerr = sum((residuals.^2));
subplot(1,2,2)
plot(tmod+mint, ymod, '-', 'LineWidth', 2)
hold on
plot(torig, ydata, 'o', 'LineWidth', 2)
legend('normal model fit', 'data')
legend boxoff
set(gca,'FontSize',16,'LineWidth',1.5)
xlabel ('days since symptom onset')
ylabel('Density')
title(['Flu normal fit sum-squared error =', num2str(sumsqerr)])

save('../out/flu_normparams.mat', 'pbestnorm')



%% Fit SARS to normal and gamma distributions
disease = 'SARS';
distrib = 'gamma';
ydata = figdigHeSARS(:,2);
ydata = ydata./sum(ydata);
%ydata = ydata- min(ydata);
torig = figdigHeSARS(:,1);
t = torig-min(torig);
p0 = [10, 0.7];
yguess = gampdf(torig, p0(1), p0(2)); 

% Fit SARS to gamma 
err = @(p) objfun(p,t, ydata,disease, distrib);
[pbest, resnorm, residuals] = lsqnonlin(err, p0);
tmod = t(1):0.1:t(end);
ymod = gampdf(tmod, pbest(1), pbest(2));
sumsqerr = sum((residuals.^2));
figure;
subplot(1,2,1)
plot(tmod+min(torig), ymod, '-', 'LineWidth', 2)
hold on
%plot(torig, yguess, '.', 'LineWidth',2)
plot(torig, ydata, 'o', 'LineWidth', 2)
legend('gamma model fit', 'data')
legend boxoff
set(gca,'FontSize',16,'LineWidth',1.5)
xlabel ('days since symptom onset')
ylabel('Density')
title(['SARS gamma fit sum-squared error =', num2str(sumsqerr)])
save('../out/SARS_gammaparams.mat', 'pbest')

disease = 'SARS';
distrib = 'normal';
ydata = figdigHeSARS(:,2);
ydata = ydata./sum(ydata);
ydata = ydata- min(ydata);
torig = figdigHeSARS(:,1);
p0 = [10, 2];
yguess = normpdf(torig, p0(1), p0(2)); 

% Fit flu to normal
err = @(p) objfun(p,torig, ydata,disease, distrib);
[pbestnorm, resnorm, residuals] = lsqnonlin(err, p0);
tmod = torig(1):0.1:torig(end);
ymod = normpdf(tmod, pbestnorm(1), pbestnorm(2));
sumsqerr = sum((residuals.^2));

subplot(1,2,2)
plot(tmod, ymod, '-', 'LineWidth', 2)
hold on
%plot(torig, yguess, '.', 'LineWidth', 2)
plot(torig, ydata, 'o', 'LineWidth', 2)
legend('normal model fit','data')
legend boxoff
set(gca,'FontSize',16,'LineWidth',1.5)
xlabel ('days since symptom onset')
ylabel('Density')
title(['SARS normal fit sum-squared error =', num2str(sumsqerr)])

save('../out/SARS_normparams.mat', 'pbestnorm')




%% Fit flu to normal
distrib = 'normal';
err = @(p) objfun(p,t, ydata,disease, distrib, weights);
[pbest, resnorm, residuals] = lsqnonlin(err, p0);
ymod = normpdf(t, pbest(1), pbest(2));
sumsqerr = sum((residuals.^2));
figure;
plot(torig, ymod, '*')
hold on
plot(torig, ydata, 'o')
legend('normal model fit', 'data')
set(gca,'FontSize',16,'LineWidth',1.5)
xlabel ('days since symptom onset')
ylabel('Density')
title(['sum-squared error =', num2str(sumsqerr)])


