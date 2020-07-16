%SEI-Isr-R Model for FLU!! 
% to evaluate effect of test-isolation strategy
% Idea is to use the proposed model structure and to tweak the
% parameters/infectiousness function to determine the region in which a
% test-isolation strategy would be reasonably effective for SARS-CoV-2 (and
% compare it to SARS-CoV-1):
% From Ip et al analysis of flu from 2008- 2014 we know:
% infectiousness peaks right at symptom onset
% and that there are ~11% asymptomatic cases, 13% pauci symptomatic, but
% mean viral levels were 1-2log10 levels lower than symptomatic
% typically infectious 2 days before and peak one day after symptom onset
% and continue to be infectious until 6 days after symptom onset

% ideal isolation: 1 day after symptom onset (3 days after infection), 76% isolated
% realistic isolation: 2 days after symptom onset ( 4.5 days after infection), 50% isolated
 


% Model structure:
% dS/dt = -beta(tau)*S*I;
% dE/dt = beta(tau)*S*I-gamma*E;
% dI/dt = gamma*E - alpha*delta*I -sigma*I
% dIsr/dt = alpha*delta*I -;
% dR/dt = sigma*I;

%June 19 2020
close all; clear all; clc;
%% Make plot of time-dependent infectivity 
% pick some gamma hyperparameters that match figure from Ip et al
dtau = 0.1;
tend = 200;
tvec = 0:dtau:tend;
tvec = tvec';
tsym = 1.87;
a = 3.87; % 
b = 0.945; % Guess these to match Ip et al figure
beta0 = 2.3; 
w_tau = gampdf(tvec, a,b);
disease = 'flu';

figure;
plot(tvec-tsym, beta0.*w_tau, '-', 'LineWidth', 2)
hold on
plot([0, 0], [0 0.28*beta0], 'k--')
legend('\beta(\tau)', 'day of symptom onset')
legend('boxoff')
ylim([0 0.28*beta0])
xlim([ -2 8])
set(gca,'FontSize',16,'LineWidth',1.5)
xlabel ('time since infection (days)')
ylabel ('Infectiousness')
title('Flu')



AUC=sum(w_tau)*dtau*beta0
%% Run SEIR model with & without isolation 
% Set IC to a single exposure
S0 = 999;
E0 = 0;
I0 = 1;
Isr0 = 0;
R0 = 0;

y0 = [S0, E0, I0, Isr0, R0];
N = sum(y0);
beta0 = 2.5; % R0 for 1918 flu
gamma = 1/1.6; % 1/gamma = duration exposed but not detectable or infectious
alpha = 0; % alpha = proportion who will become symptomatic and seek testing, 0.76 for ideal
infend = 10; % time to be considered "recovered"
trem = 1+tsym; % time from infeciton to be removed
params = [beta0, gamma, alpha, infend, trem];
tau_params = horzcat(a,b);
dt = 0.1;
tvec = 0:dt:tend;

%% Run ideal isolated and non-isolated scenarios
[y, B, new_inf, beta_t, inf_distrib, Reffn] = fwd_SEIRD_model(params,tau_params, tvec, y0, dt, disease);

ihalf = find(y(:,1)<0.8*N, 1, 'first');
thalf = tvec(ihalf);
totinf = N-y(end,1);
pctinf = 100*totinf/N;
inf = y(:,3)+y(:,4);
peakinf = 100*max(inf)/N;

colorsets1 = varycolor(500);
colorsets = varycolor(7);
figure;
for i = 1:500
    
plot(tvec, inf_distrib(:,i), '-','color', colorsets1(i,:), 'LineWidth', 2)
hold on
AUCi(round(i/100+1,0))=sum(inf_distrib(:,i)); 
%text(tauvec(i), inf_distrib(i, i), ['day ', num2str(i*dt)])
end
%legend('day 0', 'day 1', 'day 2', 'day 3', 'day 4', 'day 5', 'day 6')
%legend boxoff
xlim([ 0 14])
xlabel('infection time (\tau) (days)')
ylabel('Number of individuals X infectiousness')
title('Weighted infectiousness if no isolation')
set(gca,'FontSize',16,'LineWidth',1.5)
inf_distribNI=inf_distrib(:,500); % save for comparison figure
inf_distrib1 = inf_distrib;

figure;
subplot(2,2,1)
for i = 1:length(y0)
    plot(tvec, y(:,i), '-', 'LineWidth',2)
    hold on
end
%plot(thalf, 0.8*N, 'k*', 'LineWidth', 4)
%text(thalf+5, 0.8*N, [num2str(round(thalf,0)), ' days to 20% infected'], 'FontSize', 14)
%text(tvec(end-150),y(end-150,1), [num2str(round(pctinf,0)), '% infected'], 'FontSize',14)
%ylim([ 0 0.1*N])
legend( 'S','E', 'I', 'I_{sr}', 'R', 't_{20% infected}', 'Location', 'NorthWest')
legend boxoff
xlabel('epidemic time (days)')
ylabel('Number of individuals')
set(gca,'FontSize',16,'LineWidth',1.5)
%title(['No isolation'])
title(['No isolation, total infected = ', num2str(round(pctinf,0)),'%'])
%ylim([0 3.5e4])

subplot(2,2,2)
plot(tvec,new_inf/dt,'-', 'color',colorsets(6,:),'LineWidth', 2)
hold on
%plot(tvec, new_exp, 'm-', 'lineWidth', 2)
%plot(t, effI, 'LineWidth', 2)
%legend( 'new infections per day')
%legend boxoff
xlabel('epidemic time (days)')
ylabel('New infections per day')
set(gca,'FontSize',16,'LineWidth',1.5)
%title('Percent growth rate per day')
title(['No isolation ', num2str(round(peakinf,0)),'% infected at peak'])
ylim([0 500])
alpha2 = 0.76; % alpha = proportion who will become symptomatic and seek testing
infend = 10; % time to be considered "recovered"
trem = 1+tsym; % time from infeciton to be removed
params2 = [beta0, gamma, alpha2, infend, trem];
[y, B, new_inf, beta_t, inf_distrib, Reffp] = fwd_SEIRD_model(params2,tau_params, tvec, y0, dt, disease);
% ihalf = find(y(:,1)<0.8*N, 1, 'first');
% thalf = tvec(ihalf);
totinf = N-y(end,1);
pctinf = 100*totinf/N;
inf = y(:,3)+y(:,4);
peakinf = 100*max(inf)/N;
subplot(2,2,3)
for i = 1:length(y0)
    plot(tvec, y(:,i), '-', 'LineWidth',2)
    hold on
end
% plot(thalf, 0.8*N, 'k*', 'LineWidth', 4)
% text(thalf+5, 0.8*N, [num2str(round(thalf,0)), ' days to 20% infected'], 'FontSize', 14)
%text(tvec(end-150),y(end-150,1), [num2str(round(pctinf,0)), '% infected'], 'FontSize',14)
%ylim([ 0 0.1*N])
legend( 'S','E', 'I', 'I_{sr}', 'R', 't_{20% infected}', 'Location', 'NorthWest')
legend boxoff
xlabel('epidemic time (days)')
ylabel('Number of individuals')
set(gca,'FontSize',16,'LineWidth',1.5)
%title('Perfect symptomatic isolation')
title(['Perfect symptomatic isolation, total infected= ', num2str(round(pctinf,0)),'%'])
%title(['Perfect symptomatic isolation ',  num2str(round(pctinf,0)), '% infected total'])
%ylim([0 3.5e4])

subplot(2,2,4)
plot(tvec,new_inf/dt, '-', 'color', colorsets(6,:),'LineWidth', 2)
hold on
%plot(tvec, new_exp, 'm-', 'lineWidth', 2)
%plot(t, effI, 'LineWidth', 2)
%legend( 'new infections per day')
%legend boxoff
xlabel('epidemic time (days)')
ylabel('New infections per day')
set(gca,'FontSize',16,'LineWidth',1.5)
title(['Perfect symptomatic isolation ',num2str(round(peakinf,0)),'% infected at peak'])
ylim([0 500])
%% Plot of individuals at each stage of infection and 
figure;
plot(tvec, inf_distrib1(:,1), 'b-', 'LineWidth',5)  
hold on
plot(tvec, inf_distrib(:,1), 'r-', 'LineWidth', 5)
for i = 1:500
plot(tvec, inf_distrib1(:,i), 'b-', 'LineWidth',5)   
plot(tvec, inf_distrib(:,i), 'r-', 'LineWidth', 5)
hold on 
%text(tauvec(i), inf_distrib(i, i), ['day ', num2str(i*dt)])
end
%legend('day 0', 'day 1', 'day 2', 'day 3', 'day 4', 'day 5', 'day 6')
%legend boxoff
AUCNI = sum(inf_distrib1(:,500)).*dt;
AUCPSI = sum(inf_distrib(:,500)).*dt;
legend(['no isolation, AUC= ',num2str(round(AUCNI,0))], ['perfect symptomatic isolation, AUC= ', num2str(round(AUCPSI,1))])
legend boxoff
xlim([ 0 14])
xlabel('infection stage (\tau)(days)')
ylabel('Number of individuals x relative infectiousness')
title('Comparison of weighted infectiousness at t=50 days')
set(gca,'FontSize',16,'LineWidth',1.5)



%% Adjust one parameter at a time form the perfect symptomatic removal scenario
% Effect of increasing the time to removal
trem = tsym;
alpha = 0.76;
params = [beta0, gamma, alpha, infend, trem];
% Time to remove
figure;
tremvec = [tsym+0:2:tsym+8];
alphavec = [0:0.2:0.8];
alphavec = fliplr(alphavec);
paramsi = params;
colorsets2 = varycolor(length(tremvec));
for i = 1:length(tremvec)

    paramsi(5) = tremvec(i);
    [y, B, new_inf, R_t, new_exp,Reff] = fwd_SEIRD_model(paramsi,tau_params, tvec, y0, dt, disease);
    tot_inf = sum(y(:,3));
    S_t(:,i) = y(:,1);
    new_infi(:,i) = new_inf;
%     i20 = find(y(:,1)<0.8*N, 1, 'first');
%     t20(i) = tvec(i20);
    pct_infected(i) = 100.*(1-(S_t(end,i)./N));
    inf = y(:,3)+y(:,4);
    peakinf(i) = 100*max(inf)/N;
    Reffi(i) = Reff;
    
end

subplot(3,3,1)
for i = 1:length(tremvec)
    plot(tvec, S_t(:,i)/N, '-', 'color', colorsets2(i,:), 'LineWidth',2)
    hold on
    %text(tvec(20*i), S_t(20*i), ['t_{isolate post symptoms}=', num2str(tremvec(i)-3)], 'FontSize', 14)
    
end

legend( '0 days', '2 days', '4 days', '6 days', '8 days', 'Location', 'NorthEast') 
legend boxoff
xlabel('epidemic time(days)')
ylabel('Proportion not infected')
set(gca,'FontSize',14,'LineWidth',1.5)
%title(' infections over time')
ylim([0 1])

subplot(3,3,2)
for i = 1:length(tremvec)
    plot(tvec, new_infi(:,i)./dt, '-', 'color', colorsets2(i,:), 'LineWidth',2)
    hold on
    %text(tvec(20*i), 100*new_infi(20*i, i)/N, ['t_{isolate post symptoms}=', num2str(tremvec(i)-3)], 'FontSize', 14)

end
legend( '0 days', '2 days', '4 days', '6 days', '8 days') 
legend boxoff
xlabel('epidemic time(days)')
ylabel('New infections per day')
set(gca,'FontSize',14,'LineWidth',1.5)
%title('Effect of removal time on growth rate')
ylim([ 0 500])

subplot(3,3,3)
for i = 1:length(tremvec)
    plot(tremvec(i)-tsym, Reffi(i), '*','color', colorsets2(i,:), 'LineWidth', 10)
    hold on
end
%legend( '0 days', '1 day', '2 days', '3 days', '4 days', '5 days','Location', 'NorthWest')
%legend boxoff
ylim([0 3])
xlabel('Days after symptom onset an individual is isolated')
ylabel('R_{eff}')
set(gca,'FontSize',14,'LineWidth',1.5)
%title('Effect of removal time peak infections')


% Effect of reducing the % removed
trem = tsym;
alpha = 0.76;
params = [beta0, gamma, alpha, infend, trem];
% Time to remove

paramsi = params;

for i = 1:length(alphavec)

    paramsi(3) = alphavec(i);
    [y, B, new_inf, R_t, new_exp, Reff] = fwd_SEIRD_model(paramsi,tau_params, tvec, y0, dt, disease);
    tot_inf = sum(y(:,3));
    S_t(:,i) = y(:,1);
    new_infi(:,i) = new_inf;
    pct_infected(i) = 100.*(1-(S_t(end,i)./N));
    inf = y(:,3)+y(:,4);
    peakinf(i) = 100*max(inf)/N;
    Reffi(i) = Reff;
    
end

subplot(3,3,4)
for i = 1:length(alphavec)
    plot(tvec, S_t(:,i)/N, '-', 'color', colorsets2(i,:), 'LineWidth',2)
    hold on
    %text(tvec(20*i), S_t(20*i), ['t_{isolate post symptoms}=', num2str(tremvec(i)-3)], 'FontSize', 14)
    
end
legend(  '80% ', '60%', '40%', '20%', '0%', 'Location', 'NorthWest')
legend boxoff
xlabel('epidemic time(days)')
ylabel('Proportion not infected')
set(gca,'FontSize',14,'LineWidth',1.5)
%title(' infections over time')
ylim([0 1])

subplot(3,3,5)
for i = 1:length(alphavec)
    plot(tvec, new_infi(:,i)./dt, '-', 'color', colorsets2(i,:), 'LineWidth',2)
    hold on
    %text(tvec(20*i), 100*new_infi(20*i, i)/N, ['t_{isolate post symptoms}=', num2str(tremvec(i)-3)], 'FontSize', 14)

end
legend(  '80% ', '60%', '40%', '20%', '0%')
legend boxoff
xlabel('epidemic time(days)')
ylabel('New infections per day')
set(gca,'FontSize',14,'LineWidth',1.5)
ylim([ 0 500])
%title('Effect of removal time on growth rate')

subplot(3,3,6)
for i = 1:length(alphavec)
    plot(100*alphavec(i), Reffi(i), '*','color', colorsets2(i,:), 'LineWidth', 10)
    hold on
end
ylim([0 3])
xlabel('Percent of cases isolated')
ylabel('R_{eff}')
set(gca,'FontSize',14,'LineWidth',1.5)


%Ideal vs. Realistic scenario
tremvec = [3, 4, 3];
alphavec = [0.76, 0.5, 0];
colorsets3 = [0 0 1;0 1 0; 1 0 0];
paramsi = params;
for i = 1:length(alphavec)
    paramsi(3) = alphavec(i);
    paramsi(5) = tremvec(i);
    [y, B, new_inf, R_t, new_exp, Reff] = fwd_SEIRD_model(paramsi,tau_params, tvec, y0, dt, disease);
    tot_inf = sum(y(:,3));
    S_t(:,i) = y(:,1);
    new_infi(:,i) = new_inf;
    %ihalf = find(y(:,1)<0.5*N, 1, 'first');
    %thalf(i) = tvec(ihalf);
    pct_infected(i) = 100.*(1-(S_t(end,i)./N));
    inf = y(:,3)+y(:,4);
    peakinf(i) = 100*max(inf)/N;
    Reffi(i) = Reff;
    
end
% Ideal vs realistic

subplot(3,3,7)
for i = 1:length(alphavec)
plot(tvec, S_t(:,i)/N, '-','color', colorsets3(i,:), 'LineWidth', 2)
hold on
end
legend('perfect', 'realistic','no isolation')
legend boxoff
ylabel('Percent not infected')
xlabel('epidemic time (days)')
set(gca,'FontSize',14,'LineWidth',1.5)
ylim([0 1])
%%title(['Ideal ', num2str(round(pct_infected(1),0)),'% infected vs realistic ', num2str(round(pct_infected(2),0)),'% infected scenarios'])

subplot(3,3,8)
for i=1:length(alphavec)
plot(tvec, new_infi(:,i)/dt, '-', 'color', colorsets3(i,:), 'LineWidth',2)
hold on
end
legend('perfect', 'realistic','no isolation')
legend boxoff
ylabel('New infections per day')
xlabel('epidemic time (days)')
set(gca,'FontSize',14,'LineWidth',1.5)
ylim([ 0 500])
%title('Ideal vs. realistic vs. no isolation epidemic growth rate dynamics')
scenarios = {'perfect symptomatic isolation realistic no isolation'};
subplot(3,3,9)
for i = 1:length(alphavec)
bar( i, Reffi(i), 'FaceColor', colorsets3(i,:))
hold on
end
%legend('perfect', 'realistic','no isolation')
%legend boxoff
set(gca,'FontSize',14,'LineWidth',1.5, 'XTickLabel', [])
ylabel ('R_{eff}')
ylim([0 3])
%% Heatmap of % removed vs time of removal colored by total % infected
alphav = [0:.1:0.9]; % vary % that gets isolated
alphavdense = [0:0.01:0.9];
tremv = [0:1:9]; % vary the time of removal
tremvdense = [0:0.1:9];
[ALPHA,TREM] = meshgrid(alphav,tremv); % big ol grid of parameters
ALPHAflat = reshape(ALPHA,1,[]);
TREMflat = reshape(TREM,1, []);
paramsi = params;
y0i = y0;
I0vec =[0.1 1 10 100];
sdvec = [0.5 0.75 1 1.5];
beta0vec = sdvec.*beta0;
cell{1}.I0 = y0(3);
%% Run a loop to save the time to reach 20% infections distribuion and the
% total number infected distribution
for i = 1:length(beta0vec)
for j = 1:length(ALPHAflat)
    paramsi(1) = beta0vec(i);
    paramsi(3) = ALPHAflat(j);
    paramsi(5) = TREMflat(j);
    [y, B, new_inf, R_t, new_exp, Reff] = fwd_SEIRD_model(paramsi,tau_params, tvec, y0i, dt, disease);
    tot_inf(j) = sum(y(:,3));
    new_infi(:,j) = new_inf;
%     iquart = find(y(:,1)<0.8*N, 1, 'first');
%     if isnan(iquart)
%         tquart(j) = 0;
%     end
%     if ~isnan(iquart)
%         tquart(j) = tvec(iquart);
%     end
    pct_infected(j) = 100.*(1-(y(end,1)./N));
     inf = y(:,3)+y(:,4);
    peakinf(j) = 100*max(inf)/N;
    Reffj(j) = Reff;

end
 PCTINF = reshape(pct_infected, size(ALPHA));
 %TQUART = reshape(tquart, size(ALPHA));
 PEAKINF = reshape(peakinf, size(ALPHA));
 REFF = reshape(Reffj, size(ALPHA));
 cell{i}.I0 = y0(3);
 cell{i}.R0 = beta0vec(i);
 cell{i}.PCTINF = PCTINF;
 %cell{i}.TQUART = TQUART;
 cell{i}.PEAKINF = PEAKINF;
 cell{i}.REFF = REFF;
end
%% Plot heatmaps for total infected, time to 20% infected, and peak infections 

 figure; 
for i = 1:length(beta0vec)
    
    subplot(2, 4, i)
minmax = @(x)([min(x) max(x)]);
imagesc(minmax(1*(alphav)),minmax(tremv-tsym),cell{i}.PCTINF);
[C,h]=contourf(100*ALPHA,TREM-tsym,cell{i}.PCTINF); clabel(C,h);


if i == length(beta0vec)
    colorbar;
end
%[C,h]=contourf(ALPHA,TREM,PCTINF); clabel(C,h); colorbar
hold on
%plot(100*alphav, 0*ones(length(alphav),1), 'w--', 'LineWidth', 1.5)
%plot([76 76], [-tsym 9-tsym], 'w--', 'LineWidth', 1.5)
plot([0 76], [0 0], 'g-', 'LineWidth', 3)
plot([76 76], [0 9-tsym], 'g-', 'LineWidth', 3)
plot([0 76], [9-tsym 9-tsym], 'g-', 'LineWidth', 3)
plot( [0 0], [0 9-tsym], 'g-', 'LineWidth', 3)
caxis([0 100]);


colormap(jet); 

xlabel('% of population removed');
ylabel('days after symptom onset');
%title(['R_{0}^{flu}= ',num2str(beta0vec(i)), ', total % infected']);
title([num2str(sdvec(i)*100), '% transmission, Total % Infected'])
%title([ num2str(sdvec(i)), 'R_{0}^{flu}, total % infected']);
set(gca,'Ydir','normal'); hold on;
set(gca,'FontSize',16,'LineWidth',1.5)
end

for i =1:length(beta0vec)
    subplot(2,4,i+4)
minmax = @(x)([min(x) max(x)]);
imagesc(minmax(1*(alphav)),minmax(tremv-tsym),cell{i}.REFF);
[C,h]=contourf(100*ALPHA,TREM-tsym,cell{i}.REFF); clabel(C,h); 
hold on
[x,y,z] = C2xyz(C);
j=find(z==1);
plot(x{j}, y{j}, 'w-', 'LineWidth', 3)
if i == length(beta0vec)
    colorbar;
end
%[C,h]=contourf(ALPHA,TREM,PCTINF); clabel(C,h); colorbar
hold on
%plot(100*alphav, 0*ones(length(alphav),1), 'w--', 'LineWidth', 1.5)
%plot([76 76], [-tsym 9-tsym], 'w--', 'LineWidth', 1.5)
plot([0 76], [0 0], 'g-', 'LineWidth', 3)
plot([76 76], [0 9-tsym], 'g-', 'LineWidth', 3)
plot([0 76], [9-tsym 9-tsym], 'g-', 'LineWidth', 3)
plot( [0 0], [0 9-tsym], 'g-', 'LineWidth', 3)
caxis([0 4]);

colormap(jet); 

xlabel('% of population removed');
ylabel('days after symptom onset');
title([num2str(sdvec(i)*100), '% transmission, R_{eff}'])
%title(' Peak % infected');
set(gca,'Ydir','normal'); hold on;
set(gca,'FontSize',16,'LineWidth',1.5)
end