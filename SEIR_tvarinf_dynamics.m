%SEI-Isr-R Model to evaluate effect of test-isolation strategy
% Idea is to use the proposed model structure and to tweak the
% parameters/infectiousness function to determine the region in which a
% test-isolation strategy would be reasonably effective for SARS-CoV-2 (and
% compare it to SARS-CoV-1

% Model structure:
% dS/dt = -beta(tau)*S*I;
% dE/dt = beta(tau)*S*I-gamma*E;
% dI/dt = gamma*E - alpha*delta*I -sigma*I
% dIsr/dt = alpha*delta*I;
% dR/dt = sigma*I;

%June 8 2020
close all; clear all; clc;
%% Make plot of time-dependent infectivity 
dt = 0.1;
tend = 200;
tvec = 0:dt:tend;
tvec = tvec';
a = 2.12; % from He et al
b = 1/0.69; % from He et al 
beta0 = 2.5; 
w_tau = gampdf(tvec, a,b);
tsym = 2.3; 
figure;
plot([0, 0], [0 0.28*beta0], 'k--')
hold on
plot(tvec-tsym, beta0.*w_tau, '-', 'LineWidth', 2)
legend('symptom onset')
legend('boxoff')
ylim([0 0.28*beta0])
xlim([ -tsym 11])
set(gca,'FontSize',16,'LineWidth',1.5)
xlabel ('Days after symptom onset')
ylabel ('Infectiousness')
title('SARS-CoV-2')


D = sum(w_tau.*dt.*tvec)

AUC=sum(w_tau)*dt*beta0
%% Run SEIR model with & without isolation 
% Set IC to a single exposure
S0 = 999;
E0 = 0;
I0 = 1;
Isr0 = 0;
R0 = 0;

y0 = [S0, E0, I0, Isr0, R0];
N = sum(y0);
beta0 = 2.5; % initial rate of infection
gamma = 1/3; % 1/gamma = duration exposed but not detectable or infectious
alpha = 0; % alpha = proportion who will become symptomatic and seek testing
infend = 14; % time to be considered "recovered"
trem = 0+tsym; % time from infeciton to be removed
params = [beta0, gamma, alpha, infend, trem];
tau_params = horzcat(a,b);
D = a/b;
dt = 0.1;
tvec = 0:dt:tend;
disease = 'COVID';

%% Run ideal isolated and non-isolated scenarios
[y, B, new_inf, beta_t,inf_distrib, Reff] = fwd_SEIRD_model(params,tau_params, tvec, y0, dt, disease);
ihalf = find(y(:,1)<0.8*N, 1, 'first');
thalf = tvec(ihalf);
i1001 = find(y(:,3)>30, 1, 'first');
t1001 = tvec(i1001);
totinf = N-y(end,1);

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
ylabel('Force of infection (I*\beta)')
title('Weighted infectiousness if no isolation')
set(gca,'FontSize',16,'LineWidth',1.5)
inf_distribNI=inf_distrib(:,500); % save for comparison figure
inf_distrib1 = inf_distrib;

figure;
subplot(2,2,1)
for i = 1:length(y0)
    plot(tvec, 100*y(:,i)./N, '-', 'LineWidth',2)
    %plot(tvec, 100*y(:,i)./N, '-', 'color', colorsets(i+5,:), 'LineWidth',2)
    hold on
end
%plot(thalf, 80, 'k*', 'LineWidth', 4)
%text(thalf+5, 80, [num2str(round(thalf,0)), ' days to 20% infected'], 'FontSize', 14)
%text(tvec(end-150),y(end-150,1), [num2str(round(pctinf,0)), '% infected'], 'FontSize',14)
%ylim([ 0 0.1*N])
legend( 'S','E', 'I', 'I_{sr}', 'R', 't_{20% infected}', 'Location', 'NorthWest')
legend boxoff
xlabel('epidemic time (days)')
ylabel('% of population')
set(gca,'FontSize',16,'LineWidth',1.5)
%title(['No isolation'])
title(['No isolation, total infected = ', num2str(round(pctinf,0)),'%'])
%ylim([0 3.5e4])


subplot(2,2,2)
plot(tvec,new_inf./dt,'-', 'color',colorsets(6,:),'LineWidth', 2)
hold on
%plot(tvec, new_exp, 'm-', 'lineWidth', 2)
%plot(t, effI, 'LineWidth', 2)
%legend( 'new infections per day')
%legend boxoff
xlabel('epidemic time (days)')
ylabel('New infections per day')
set(gca,'FontSize',16,'LineWidth',1.5)
%title('Percent growth rate per day')
ylim([ 0 400])
title(['No isolation ', num2str(round(peakinf,0)),'% infected at peak'])


% subplot(2,2,2)
% plot(tvec,Rt1,'-', 'color',colorsets(6,:),'LineWidth', 2)
% hold on
% plot(tvec, beta0*ones(length(tvec),1), 'k--', 'LineWidth', 2)
% %plot(tvec, new_exp, 'm-', 'lineWidth', 2)
% %plot(t, effI, 'LineWidth', 2)
% %legend( 'new infections per day')
% %legend boxoff
% %xlim([ 0 50])
% xlabel('epidemic time (days)')
% ylabel('R_t')
% set(gca,'FontSize',16,'LineWidth',1.5)
% %title('Percent growth rate per day')
% title(['No isolation ', num2str(round(peakinf,0)),'% infected at peak'])
%ylim([0 4e4])
alpha2 = 0.8; % alpha = proportion who will become symptomatic and seek testing
infend = 14; % time to be considered "recovered"
trem = tsym; % time from infeciton to be removed
params2 = [beta0, gamma, alpha2, infend, trem];
[y, B, new_inf, beta_t,inf_distrib, Reffp] = fwd_SEIRD_model(params2,tau_params, tvec, y0, dt,disease);
ihalf = find(y(:,1)<0.8*N, 1, 'first');
thalf = tvec(ihalf);
i100 = find(y(:,3)>30, 1, 'first');
t100 = tvec(i100);
totinf = N-y(end,1);
pctinf = 100*totinf/N;
inf = y(:,3)+y(:,4);
peakinf = 100*max(inf)/N;


subplot(2,2,3)
for i = 1:length(y0)
    plot(tvec, 100*y(:,i)./N, '-', 'LineWidth',2)
    %plot(tvec, 100*y(:,i)./N, '-', 'color', colorsets(i,:), 'LineWidth',2)
    hold on
end
%plot(thalf, 80, 'k*', 'LineWidth', 4)
%text(thalf+5, 80, [num2str(round(thalf,0)), ' days to 20% infected'], 'FontSize', 14)
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
ylim ([0 400])
ylabel('New infections per day')
set(gca,'FontSize',16,'LineWidth',1.5)
title(['Perfect symptomatic isolation ',num2str(round(peakinf,0)),'% infected at peak'])

% subplot(2,2,4)
% plot(tvec,Rt2,'-', 'color',colorsets(6,:),'LineWidth', 2)
% hold on
% plot(tvec, beta0*ones(length(tvec),1), 'k--', 'LineWidth', 2)
% %plot(tvec, new_exp, 'm-', 'lineWidth', 2)
% %plot(t, effI, 'LineWidth', 2)
% %legend( 'new infections per day')
% %legend boxoff
% xlabel('epidemic time (days)')
% ylabel('R_t')
% %xlim ([0 50])
% set(gca,'FontSize',16,'LineWidth',1.5)
% title(['Perfect symptomatic isolation ',num2str(round(peakinf,0)),'% infected at peak'])



%% Plot of individuals at each stage of infection and 
figure;
plot(tvec, inf_distrib1(:,i1001), 'b-', 'LineWidth',5)  
hold on
plot(tvec, inf_distrib(:,i100), 'r-', 'LineWidth', 5)
% for i = 1:500
% plot(tvec, inf_distrib1(:,i), 'b-', 'LineWidth',5)   
% plot(tvec, inf_distrib(:,i), 'r-', 'LineWidth', 5)
% hold on 
% %text(tauvec(i), inf_distrib(i, i), ['day ', num2str(i*dt)])
% end
%legend('day 0', 'day 1', 'day 2', 'day 3', 'day 4', 'day 5', 'day 6')
%legend boxoff
AUCNI = sum(inf_distrib1(:,500)).*dt;
AUCPSI = sum(inf_distrib(:,500)).*dt;
legend( ['no isolation, t=', num2str(t1001),' days'], ['perfect symptomatic isolation, t=', num2str(t100),' days'])
%legend(['no isolation, AUC= ',num2str(round(AUCNI,0))], ['perfect symptomatic isolation, AUC= ', num2str(round(AUCPSI,0))])
legend boxoff
xlim([ 0 14])
xlabel('infection stage (\tau)(days)')
ylabel('Number of individuals x relative infectiousness')
title('Comparison of weighted infectiousness at 30 infections')
set(gca,'FontSize',16,'LineWidth',1.5)



%% Adjust one parameter at a time form the perfect symptomatic removal scenario
% Effect of increasing the time to removal
trem = tsym+0;
alpha = 0.8;
params = [beta0, gamma, alpha, infend, trem];
% Time to remove
figure;
tremvec = [tsym:1:4+tsym];
alphavec = [0:0.2:0.8];
alphavec = fliplr(alphavec);
paramsi = params;
colorsets2 = varycolor(length(tremvec));
for i = 1:length(tremvec)

    paramsi(5) = tremvec(i);
    [y, B, new_inf, R_t, new_exp, Reff] = fwd_SEIRD_model(paramsi,tau_params, tvec, y0, dt, disease);
    tot_inf = sum(y(:,3));
    S_t(:,i) = y(:,1);
    new_infi(:,i) = new_inf;
    i20 = find(y(:,1)<0.8*N, 1, 'first');
    t20(i) = tvec(i20);
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
legend( '0 days',  '1 day',  '2 days',  '3 days', '4 days') 
legend boxoff

xlabel('epidemic time (days)')
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
legend( '0 days',  '1 day',  '2 days',  '3 days', '4 days') 
legend boxoff
xlabel('epidemic time (days)')
ylabel('New infections per day')
set(gca,'FontSize',14,'LineWidth',1.5)
%title('Effect of removal time on growth rate')

subplot(3,3,3)
for i = 1:length(tremvec)
    plot(tremvec(i)-tsym, Reffi(i), '*','color', colorsets2(i,:), 'LineWidth', 10)
    hold on
end
%legend( '0 days', '1 day', '2 days', '3 days', '4 days', '5 days','Location', 'NorthWest')
%legend boxoff
%ylim([15 60])
ylim([1 3])
xlabel('Days after symptom onset an individual is isolated')
ylabel('R_{eff}')
set(gca,'FontSize',14,'LineWidth',1.5)

% subplot(3,3,3)
% for i = 1:length(tremvec)
%     plot(tremvec(i)-tsym, peakinf(i), '*','color', colorsets2(i,:), 'LineWidth', 10)
%     hold on
% end
% %legend( '0 days', '1 day', '2 days', '3 days', '4 days', '5 days','Location', 'NorthWest')
% %legend boxoff
% %ylim([15 60])
% xlabel('Days after symptom onset an individual is isolated')
% ylabel('Peak % infected at any time')
% set(gca,'FontSize',14,'LineWidth',1.5)
% %title('Effect of removal time peak infections')


% Effect of reducing the % removed
trem = tsym;
alpha = 0.8;
params = [beta0, gamma, alpha, infend, trem];
% Time to remove

paramsi = params;

for i = 1:length(alphavec)

    paramsi(3) = alphavec(i);
    [y, B, new_inf, R_t, new_exp, Reff] = fwd_SEIRD_model(paramsi,tau_params, tvec, y0, dt, disease);
    tot_inf = sum(y(:,3));
    S_t(:,i) = y(:,1);
    new_infi(:,i) = new_inf;
    i20 = find(y(:,1)<0.8*N, 1, 'first');
    t20(i) = tvec(i20);
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
legend( '80%', '60%', '40%', '20%',  '0%')
legend boxoff

xlabel('epidemic time (days)')
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
legend( '80%', '60%', '40%', '20%',  '0%')
legend boxoff
xlabel('epidemic time (days)')
ylabel('New infections per day')
set(gca,'FontSize',14,'LineWidth',1.5)
%title('Effect of removal time on growth rate')
subplot(3,3,6)
for i = 1:length(alphavec)
    plot(100*alphavec(i), Reffi(i), '*','color', colorsets2(i,:), 'LineWidth', 10)
    hold on
end
%ylim([15 60])
ylim([1 3])
xlabel('Percent of cases isolated')
ylabel('R_{eff}')
set(gca,'FontSize',14,'LineWidth',1.5)

% subplot(3,3,6)
% for i = 1:length(alphavec)
%     plot(100*alphavec(i), peakinf(i), '*','color', colorsets2(i,:), 'LineWidth', 10)
%     hold on
% end
% %ylim([15 60])
% xlabel('Percent of cases isolated')
% ylabel('Peak % infected at any time')
% set(gca,'FontSize',14,'LineWidth',1.5)


%Ideal vs. Realistic scenario
tremvec = [tsym, 1.5+tsym, tsym];
alphavec = [0.8, 0.75*0.8, 0]; % assume you can isolate 75% of symptomatics
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
    Reffi(i) = Reff
    
end
% Ideal vs realistic

subplot(3,3,7)
for i = 1:length(alphavec)
plot(tvec, S_t(:,i)/N, '-','color', colorsets3(i,:), 'LineWidth', 2)
hold on
end
%legend('perfect', 'realistic','no isolation')
%legend boxoff
ylabel('Percent not infected')
xlabel('epidemic time (days)')
set(gca,'FontSize',14,'LineWidth',1.5)
%%title(['Ideal ', num2str(round(pct_infected(1),0)),'% infected vs realistic ', num2str(round(pct_infected(2),0)),'% infected scenarios'])

subplot(3,3,8)
for i=1:length(alphavec)
plot(tvec, new_infi(:,i)./dt, '-', 'color', colorsets3(i,:), 'LineWidth',2)
hold on
end
legend('perfect', 'realistic','no isolation')
legend boxoff
ylabel('New infections per day')
xlabel('epidemic time (days)')
set(gca,'FontSize',14,'LineWidth',1.5)
%title('Ideal vs. realistic vs. no isolation epidemic growth rate dynamics')
scenarios = {'perfect symptomatic isolation realistic no isolation'};

subplot(3,3,9)
for i = 1:length(alphavec)
bar( i, Reffi(i), 'FaceColor', colorsets3(i,:))
hold on
end
ylim([1 3])
%legend('perfect', 'realistic','no isolation')
%legend boxoff
set(gca,'FontSize',14,'LineWidth',1.5, 'XTickLabel', [])
ylabel ('R_{eff}')
%ylim([20 60])

% subplot(3,3,9)
% for i = 1:length(alphavec)
% bar( i, peakinf(i), 'FaceColor', colorsets3(i,:))
% hold on
% end
% %legend('perfect', 'realistic','no isolation')
% %legend boxoff
% set(gca,'FontSize',14,'LineWidth',1.5, 'XTickLabel', [])
% ylabel ('Peak % infected')
% %ylim([20 60])


%% What is the critical R0 to make realistic isolation reasonably effective?
% perfect, realistic, none
tremvec = [tsym, 1.5+tsym, tsym];
alphavec = [0.8, 0.75*0.8, 0]; % assume you can isolate 75% of symptomatics
colorsets3 = [0 0 1;0 1 0; 1 0 0];
R0vec = [0.5 0.625 0.75 0.875 1].*beta0;
R0vec = linspace(0, 1.2, 100).*beta0;
paramsi = params;
for j = 1:length(R0vec)
    paramsi(1) = R0vec(j);
for i = 1:length(alphavec)
    paramsi(3) = alphavec(i);
    paramsi(5) = tremvec(i);
    [y, B, new_inf, R_t, new_exp, Reff] = fwd_SEIRD_model(paramsi,tau_params, tvec, y0, dt, disease);
    tot_inf = sum(y(:,3));
    S_t(:,i) = y(:,1);
    new_infi(:,i) = new_inf;
    %ihalf = find(y(:,1)<0.5*N, 1, 'first');
    %thalf(i) = tvec(ihalf);
    pct_infected(i,j) = 100.*(1-(S_t(end,i)./N));
    inf = y(:,3)+y(:,4);
    peakinf(i,j) = 100*max(inf)/N;
    Reffi(i,j) = Reff;
    
end

end

%%
pct_trans = linspace(0.5, 1.2, 100);
ind = find(pct_infected(2,:)>=20,1, 'first');
Rcrit = R0vec(ind);
pct_red_crit = pct_trans(ind);
Reffcrit = Reffi(2, ind);

i2 = find(pct_infected(1,:)>=20, 1, 'first');
Rcritperf = R0vec(i2);
pct_red_perf_crit = pct_trans(i2);
Reffperf = Reffi(1,i2);

i5 = find(pct_infected(3,:)>=20, 1, 'first');
Rcritno = R0vec(i5);
pct_red_no_crit = pct_trans(i5);
Reffno = Reffi(1,i5);

i3 = find(Reffi(2,:)>=1, 1, 'first');
Rcrit3 = R0vec(i3);
i4 = find(Reffi(1,:)>=1, 1, 'first');
Rcrit4 = R0vec(i4);
figure;
for i = 1:3
plot(R0vec, pct_infected(i,:), 'color', colorsets3(i,:), 'LineWidth', 2)
hold on
end
hold on
plot( R0vec, 20*ones(length(R0vec),1), 'k--', 'LineWidth', 2)
%plot([Rcrit Rcrit], [0 20], 'k--', 'LineWidth', 2)
text(Rcrit, 23, ['Realistic R_{crit} =' num2str(round(Rcrit,1))], 'FontSize', 12)
%plot([Rcritperf Rcritperf], [0 20], 'k--', 'LineWidth',2)
text(Rcritperf, 23, ['Perfect R_{crit} =' num2str(round(Rcritperf,1))], 'FontSize', 12)
plot([beta0 beta0], [0 20], 'k--', 'LineWidth', 2)
text(beta0, 23, ['COVID-19 R_0 =' num2str(round(beta0,1))], 'FontSize', 12)
%plot([Rcritno Rcritno], [0 20], 'k--', 'LineWidth',2)
text(Rcritno-0.1, 17, ['No isolation R_0 =' num2str(round(Rcritno,1))], 'FontSize', 12)
legend ('perfect', 'realistic', 'no isolation', 'Location', 'NorthWest')
legend boxoff
xlabel ('R_0')
xlim([1 R0vec(end)])
title('Critical R_0 for symptomatic isolation')
ylabel ('Total % infected')
set(gca,'FontSize',16,'LineWidth',1.5)

figure;
for i = 1:3
plot(R0vec, peakinf(i,:), 'color', colorsets3(i,:), 'LineWidth', 2)
hold on
end
hold on
%plot( R0vec, 20*ones(length(R0vec),1), 'k--', 'LineWidth', 2)
%plot([Rcrit Rcrit], [0 20], 'k--', 'LineWidth', 2)
%text(Rcrit, 23, ['Realistic R_{crit} =' num2str(round(Rcrit,1))], 'FontSize', 12)
%plot([Rcritperf Rcritperf], [0 20], 'k--', 'LineWidth',2)
%text(Rcritperf, 23, ['Perfect R_{crit} =' num2str(round(Rcritperf,1))], 'FontSize', 12)
%plot([beta0 beta0], [0 20], 'k--', 'LineWidth', 2)
%text(beta0, 23, ['COVID-19 R_0 =' num2str(round(beta0,1))], 'FontSize', 12)
legend ('perfect', 'realistic', 'no isolation', 'Location', 'NorthWest')
legend boxoff
xlabel ('R_0')
xlim([R0vec(1) R0vec(end)])
title('Peak % infected vs. R_0')
ylabel ('Peak % infected')
set(gca,'FontSize',16,'LineWidth',1.5)

figure;
for i = 1:3
plot(R0vec, Reffi(i,:), 'color', colorsets3(i,:), 'LineWidth', 2)
hold on
end
hold on
plot( R0vec, ones(length(R0vec),1), 'k--', 'LineWidth', 2)
plot([Rcrit3 Rcrit3], [0 1], 'k--', 'LineWidth', 2)
text(Rcrit3, 0.9, ['Realistic R_{crit} =' num2str(round(Rcrit3,1))], 'FontSize', 12)
plot([Rcrit4 Rcrit4], [0 1], 'k--', 'LineWidth',2)
text(Rcrit4 , 0.9, ['Perfect R_{crit} =' num2str(round(Rcrit4,1))], 'FontSize', 12)
plot([beta0 beta0], [0 1], 'k--', 'LineWidth', 2)
text(beta0, 0.9, ['COVID-19 R_0 =' num2str(round(beta0,1))], 'FontSize', 12)
legend ('perfect', 'realistic', 'no isolation', 'Location', 'NorthWest')
legend boxoff
xlabel ('R_0')
xlim([R0vec(1) R0vec(end)])
title('Critical R_0 vs R_{eff} for symptomatic isolation')
ylabel ('R_{eff}')
set(gca,'FontSize',16,'LineWidth',1.5)

%% Plot I0 curves

trem = 2+tsym;
alpha = 0.6;
params = [beta0, gamma, alpha, infend, trem];
% Time to remove

I0vec = [0.01 0.1 1 10];
y0i = y0;

paramsi = params;
colorsets2 = varycolor(length(I0vec));

for i = 1:length(I0vec)
    y0i(3) = I0vec(i);
    y0i(1) = 1000-y0i(3);
    [y, B, new_inf, R_t, new_exp, Reff] = fwd_SEIRD_model(params,tau_params, tvec, y0i, dt, disease);
    tot_inf = sum(y(:,3));
    S_t(:,i) = y(:,1);
    new_infi(:,i) = new_inf;
    i20 = find(y(:,1)<0.8*N, 1, 'first');
    t20(i) = tvec(i20);
    pct_infected(i) = 100.*(1-(S_t(end,i)./N));
    inf = y(:,3)+y(:,4);
    peakinf(i) = 100*max(inf)/N;
    inf_t(:,i) = y(:,3);
    Reffj(i) = Reff;
   
    
end

figure;
subplot(1,2,1)
for i = 1:length(I0vec)
    plot(tvec, 100*S_t(:,i)/N, '-', 'color', colorsets(i,:), 'LineWidth',2)
    hold on
    %text(tvec(20*i), S_t(20*i), ['t_{isolate post symptoms}=', num2str(tremvec(i)-3)], 'FontSize', 14)
    
end
legend( 'I_0 = 0.001%','I_0 = 0.01%', 'I_0 = 0.1%', 'I_0 = 1%' ) 
legend boxoff

xlabel('epidemic time (days)')
ylabel('% of population not infected')
set(gca,'FontSize',14,'LineWidth',1.5)
%title(' infections over time')
subplot(1,2,2)
for i = 1:length(I0vec)
    plot(tvec, 100*inf_t(:,i)/N, '-', 'color', colorsets(i,:), 'LineWidth',2)
    hold on
    %text(tvec(20*i), S_t(20*i), ['t_{isolate post symptoms}=', num2str(tremvec(i)-3)], 'FontSize', 14)
    
end
legend( 'I_0 = 0.001%','I_0 = 0.01%', 'I_0 = 0.1%', 'I_0 = 1%' ) 
legend boxoff

xlabel('epidemic time (days)')
ylabel('% of population')
set(gca,'FontSize',14,'LineWidth',1.5)




%% Heatmap of % removed vs time of removal colored by total % infected
alphav = [0.1:.1:0.9]; % vary % that gets isolated
tremv = [0:1:8]; % vary the time of removal
[ALPHA,TREM] = meshgrid(alphav,tremv); % big ol grid of parameters
ALPHAflat = reshape(ALPHA,1,[]);
TREMflat = reshape(TREM,1, []);
paramsi = params;
y0i = y0;
I0vec =[0.1 1 10 100];
sdvec = [0.5 0.75 1 1.5];
%sdvec = [0.9 1 1.05 1.1];
beta0vec = beta0.*sdvec;
cell{1}.I0 = y0(3);
pct_infected = [];
peakinf = [];
%% Run a loop to save the total % infected and the Reff
% total number infected distribution
for i = 1:length(beta0vec)
for j = 1:length(ALPHAflat)
    paramsi(1) = beta0vec(i);
    paramsi(3) = ALPHAflat(j);
    paramsi(5) = TREMflat(j);
    [y, B, new_inf, R_t, new_exp, Reff] = fwd_SEIRD_model(paramsi,tau_params, tvec, y0i, dt,disease);
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
     Reffvec(j) = Reff;

end
 PCTINF = reshape(pct_infected, size(ALPHA));
 
 PEAKINF = reshape(peakinf, size(ALPHA));
 REFF = reshape(Reffvec, size(ALPHA));
 cell{i}.I0 = y0(3);
 cell{i}.R0 = beta0vec(i);
 cell{i}.PCTINF = PCTINF;

 cell{i}.PEAKINF = PEAKINF;
 cell{i}.REFF = REFF;
end
%% Plot heatmaps for total infected, time to 20% infected, and peak infections 

 figure; 
for i = 1:length(beta0vec)
    
    subplot(2, 4, i)
minmax = @(x)([min(x) max(x)]);
imagesc(minmax(1*(alphav)),minmax(tremv-tsym),cell{i}.PCTINF);
%imagesc(minmax(1*(alphav)),minmax(tremv),PCTINF);
[C,h]=contourf(100*ALPHA,TREM-tsym,cell{i}.PCTINF); clabel(C,h); 


if i == length(beta0vec)
    colorbar;
end
%[C,h]=contourf(ALPHA,TREM,PCTINF); clabel(C,h); colorbar
hold on
% plot(100*alphav, 0*ones(length(alphav),1), 'w--', 'LineWidth', 1.5)
% plot([80 80], [-tsym 8-tsym], 'w--', 'LineWidth', 1.5)
plot([10 80], [0 0], 'g-', 'LineWidth', 3)
plot([80 80], [0 8-tsym], 'g-', 'LineWidth', 3)
plot([10 80], [8-tsym-.05 8-tsym-.05], 'g-', 'LineWidth', 3)
plot( [11 11], [0 8-tsym-.05], 'g-', 'LineWidth', 3)
caxis([0 100]);
%redgreencmap
%cmap = [[zeros(32,1);[.2:.8/31:1]'] zeros(64,1) [1-[0:0.8/31:.8]';zeros(32,1)] ]; % green to red colormap
colormap(jet); 
xlabel('% of population removed');
ylabel('days after symptom onset');
title([num2str(sdvec(i)*100), '% transmission, Total % Infected'])
%title(['R_{0}^{COVID}= ',num2str(beta0vec(i)), ', Total % Infected']);
set(gca,'Ydir','normal'); hold on;
set(gca,'FontSize',16,'LineWidth',1.5)
end

% Reff

for i =1:length(I0vec)
    xvals = [];
    yvals = [];
    subplot(2,4,i+4)
minmax = @(x)([min(x) max(x)]);
imagesc(minmax(1*(alphav)),minmax(tremv-tsym),cell{i}.REFF);
[C,h]=contourf(100*ALPHA,TREM-tsym,cell{i}.REFF); clabel(C,h);
[x,y,z] = C2xyz(C);
j=find(z==1);
hold on
plot(x{j}, y{j}, 'w-', 'LineWidth', 3)


if i == length(beta0vec)
    colorbar;
end
hold on
%plot(100*alphav, zeros(length(tremv),1), 'w--', 'LineWidth', 1.5)
%plot([80 80], [-tsym 8-tsym], 'w--', 'LineWidth', 1.5)
plot([10 80], [0 0], 'g-', 'LineWidth', 3)
plot([80 80], [0 8-tsym], 'g-', 'LineWidth', 3)
plot([10 80], [8-tsym-.05 8-tsym-.05], 'g-', 'LineWidth', 3)
plot( [11 11], [0 8-tsym-.05], 'g-', 'LineWidth', 3)
caxis([0 4]);
%redgreencmap
%cmap = [[zeros(32,1);[.2:.8/31:1]'] zeros(64,1) [1-[0:0.8/31:.8]';zeros(32,1)] ]; % green to red colormap
colormap(jet); 
xlabel('% of population removed');
ylabel('days after symptom onset');
title([num2str(sdvec(i)*100), '% transmission, R_{eff}'])
%title(['R_{0}^{COVID}= ',num2str(beta0vec(i)), ', Peak % Infected']);
set(gca,'Ydir','normal'); hold on;
set(gca,'FontSize',16,'LineWidth',1.5)
end
%% Heatmap of days after sym onset vs lower transmission
tremv = linspace(0.5,tsym+4, 10);
R_0vec = linspace(1, 3.5,10);
R_0vec = fliplr(R_0vec);
[R0MAT,TREM] = meshgrid(R_0vec,tremv); % big ol grid of parameters
R0flat = reshape(R0MAT,1,[]);
TREMflat = reshape(TREM,1, []);
paramsi = params;
paramsi(3) = 0.6;


for j = 1:length(R0flat)
    paramsi(1) = R0flat(j);
    paramsi(5) = TREMflat(j);
    [y, B, new_inf, R_t, new_exp, Reff] = fwd_SEIRD_model(paramsi,tau_params, tvec, y0, dt,disease);
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
     Reffvec(j) = Reff;

end
 PCTINF = reshape(pct_infected, size(R0MAT));
 
 PEAKINF = reshape(peakinf, size(R0MAT));
 REFF = reshape(Reffvec, size(R0MAT));
%%
figure;
minmax = @(x)([min(x) max(x)]);
imagesc(minmax(1*(R_0vec)),minmax(tremv-tsym),PCTINF);
%imagesc(minmax(1*(alphav)),minmax(tremv),PCTINF);
[C,h]=contourf(R0MAT,TREM-tsym,PCTINF); clabel(C,h); 
caxis([0 100]);
colorbar;
colormap(jet); 
xlabel('R_0');
ylabel('days after symptom onset');
title('Total % Infected')
%title(['R_{0}^{COVID}= ',num2str(beta0vec(i)), ', Peak % Infected']);
set(gca,'Ydir','normal'); hold on;
set(gca,'FontSize',16,'LineWidth',1.5)

figure;
minmax = @(x)([min(x) max(x)]);
imagesc(minmax(1*(R_0vec)),minmax(tremv-tsym),REFF);
%imagesc(minmax(1*(alphav)),minmax(tremv),PCTINF);
[C,h]=contourf(100*R0MAT/beta0,TREM-tsym,REFF); clabel(C,h);
colorbar;
[x,y,z] = C2xyz(C);
j=find(z==1);
hold on
plot(x{j}, y{j}, 'w--', 'LineWidth', 3)
caxis([0 3.5]);
colormap(jet); 
xlabel('% transmission');
ylabel('days after symptom onset');
title('R_{eff}')
%title(['R_{0}^{COVID}= ',num2str(beta0vec(i)), ', Peak % Infected']);
set(gca,'Ydir','normal'); hold on;
set(gca,'FontSize',16,'LineWidth',1.5)

figure;
minmax = @(x)([min(x) max(x)]);
imagesc(minmax(1*(R_0vec)),minmax(tremv-tsym),REFF);
%imagesc(minmax(1*(alphav)),minmax(tremv),PCTINF);
[C,h]=contourf(R0MAT,TREM-tsym,REFF); clabel(C,h);
colorbar;
[x,y,z] = C2xyz(C);
j=find(z==1);
hold on
plot(x{j}, y{j}, 'w--', 'LineWidth', 3)
caxis([0 3.5]);
colormap(jet); 
xlabel('R_0');
ylabel('days after symptom onset');
title('R_{eff}')
%title(['R_{0}^{COVID}= ',num2str(beta0vec(i)), ', Peak % Infected']);
set(gca,'Ydir','normal'); hold on;
set(gca,'FontSize',16,'LineWidth',1.5)
%%
% PEAK INF
% for i =1:length(I0vec)
%     subplot(2,4,i+4)
% minmax = @(x)([min(x) max(x)]);
% imagesc(minmax(1*(alphav)),minmax(tremv-tsym),cell{i}.PEAKINF);
% [C,h]=contourf(100*ALPHA,TREM-tsym,cell{i}.PEAKINF); clabel(C,h); 
% if i == length(beta0vec)
%     colorbar;
% end
% hold on
% plot(100*alphav, zeros(length(tremv),1), 'w--', 'LineWidth', 1.5)
% plot([80 80], [-tsym 8-tsym], 'w--', 'LineWidth', 1.5)
% plot([10 80], [0 0], 'g-', 'LineWidth', 3)
% plot([80 80], [0 8-tsym], 'g-', 'LineWidth', 3)
% plot([10 80], [8-tsym-.05 8-tsym-.05], 'g-', 'LineWidth', 3)
% plot( [11 11], [0 8-tsym-.05], 'g-', 'LineWidth', 3)
% caxis([0 100]);
% %redgreencmap
% %cmap = [[zeros(32,1);[.2:.8/31:1]'] zeros(64,1) [1-[0:0.8/31:.8]';zeros(32,1)] ]; % green to red colormap
% colormap(jet); 
% xlabel('% of population removed');
% ylabel('days after symptom onset');
% title(' Peak % infected');
% title(['R_{0}^{COVID}= ',num2str(beta0vec(i)), ', Peak % Infected']);
% set(gca,'Ydir','normal'); hold on;
% set(gca,'FontSize',16,'LineWidth',1.5)
% end





%% Vary I0
cell2{1}.I0 = I0vec(1);

for i = 1:length(I0vec)
for j = 1:length(ALPHAflat)
    paramsi(1) = 2.3;
    paramsi(3) = ALPHAflat(j);
    paramsi(5) = TREMflat(j);
    y0i(3) = I0vec(i);
    y0i(1) = 1000-y0i(3);
    [y, B, new_inf, R_t, new_exp, Reff] = fwd_SEIRD_model(paramsi,tau_params, tvec, y0i, dt, disease);
    tot_inf(j) = sum(y(:,3));
    new_infi(:,j) = new_inf;
    Reffvec(j) = Reff;
    
    pct_infected(j) = 100.*(1-(y(end,1)./N));
     inf = y(:,3)+y(:,4);
     peakinf(j) = 100*max(inf)/N;
     Reffvec(j) = Reff;


end
 PCTINF = reshape(pct_infected, size(ALPHA));
 REFF = reshape(Reffvec, size(ALPHA));
 PEAKINF = reshape(peakinf, size(ALPHA));
 cell2{i}.I0 = I0vec(i);
 cell2{i}.R0 = 2.5;
 cell2{i}.PCTINF = PCTINF;

 cell2{i}.PEAKINF = PEAKINF;
 cell2{i}.REFF = REFF;
end

%% Plot heatmaps varying by I0
figure; 

for i = 1:length(I0vec)
    
    subplot(2, 4, i)
minmax = @(x)([min(x) max(x)]);
imagesc(minmax(1*(alphav)),minmax(tremv-tsym),cell2{i}.PCTINF);
%imagesc(minmax(1*(alphav)),minmax(tremv),PCTINF);
[C,h]=contourf(100*ALPHA,TREM-tsym,cell2{i}.PCTINF); clabel(C,h); colorbar
%[C,h]=contourf(ALPHA,TREM,PCTINF); clabel(C,h); colorbar
hold on
plot([10 80], [0 0], 'g-', 'LineWidth', 3)
plot([80 80], [0 8-tsym], 'g-', 'LineWidth', 3)
plot([10 80], [8-tsym-.05 8-tsym-.05], 'g-', 'LineWidth', 3)
plot( [11 11], [0 8-tsym-.05], 'g-', 'LineWidth', 3)
caxis([0 100]);
hold on
if i == length(beta0vec)
    colorbar;
end
hold on

%redgreencmap
%cmap = [[zeros(32,1);[.2:.8/31:1]'] zeros(64,1) [1-[0:0.8/31:.8]';zeros(32,1)] ]; % green to red colormap
colormap(jet); 
colorbar;
xlabel('% of population removed');
ylabel('days after symptom onset');
title(['I_0=', num2str(100*I0vec(i)/N), '%, Total % Infected']);
set(gca,'Ydir','normal'); hold on;
set(gca,'FontSize',16,'LineWidth',1.5)
end

% REFF

for i =1:length(I0vec)
    subplot(2,4,i+4)
minmax = @(x)([min(x) max(x)]);
imagesc(minmax(1*(alphav)),minmax(tremv-tsym),cell2{i}.REFF);
[C,h]=contourf(100*ALPHA,TREM-tsym,cell2{i}.REFF); clabel(C,h); colorbar
hold on

[x,y,z] = C2xyz(C);
j=find(z==1);
hold on
plot(x{j}, y{j}, 'w--', 'LineWidth', 3)
caxis([0 3.5]);
if i == length(beta0vec)
    colorbar;
end
hold on

plot([10 80], [0 0], 'g-', 'LineWidth', 3)
plot([80 80], [0 8-tsym], 'g-', 'LineWidth', 3)
plot([10 80], [8-tsym-.05 8-tsym-.05], 'g-', 'LineWidth', 3)
plot( [11 11], [0 8-tsym-.05], 'g-', 'LineWidth', 3)
caxis([0 3.5]);
%redgreencmap
%cmap = [[zeros(32,1);[.2:.8/31:1]'] zeros(64,1) [1-[0:0.8/31:.8]';zeros(32,1)] ]; % green to red colormap
colormap(jet); 
colorbar;
xlabel('% of population removed');
ylabel('days after symptom onset');
title(['I_0=', num2str(100*I0vec(i)/N), '%, R_{eff}']);
set(gca,'Ydir','normal'); hold on;
set(gca,'FontSize',16,'LineWidth',1.5)
end





%% Plot weighted infectiousness at different times
figure;
plot(tvec, B(1,:),'-', 'LineWidth',2)
hold on
plot(tvec, B(2/dt+1,:),'-', 'LineWidth',2)
plot(tvec, B(4/dt+1,:),'-', 'LineWidth',2)
plot(tvec, B(6/dt+1,:),'-', 'LineWidth',2)
plot(tvec,B(20/dt+1,:), 'g-', 'LineWidth',2)
plot( tvec, B(40/dt+1,:),'b-', 'LineWidth',2)
%plot(tauvec, B(10000,:), 'c-', 'LineWidth',2)
%plot(tauvec,B(20000,:), 'k-', 'LineWidth',2)
legend('day 0', 'day 2', 'day 4','day 6', 'day 20', 'day 40', 'day 100', 'day 200')
legend boxoff
xlim([ 0 14])
xlabel('infection age \tau (days)')
ylabel('density of infectious individuals')
set(gca,'FontSize',16,'LineWidth',1.5)
title('Density of infectious individuals over infection age')
%%
figure;
plot(tvec(10:end), log(2)./beta_t(10:end)./dt, 'r-', 'LineWidth',2)
xlabel('epidemic time (days)')
ylabel('Epidemic Doubling Time')
set(gca,'FontSize',16,'LineWidth',1.5)

figure;
plot(tvec, y(:,1)/N, 'b-', 'LineWidth', 2)
xlabel('epidemic time (days)')
ylabel('Percent not infected')
set(gca,'FontSize',16,'LineWidth',1.5)

colorsets = varycolor(length(tvec));
figure;
for i = 1:length(tvec)
    
plot(tvec, inf_distrib(:,i), '-','color', colorsets(i,:), 'LineWidth', 2)
hold on
AUCi(round(i/100+1,0))=sum(inf_distrib(:,i)); 
%text(tauvec(i), inf_distrib(i, i), ['day ', num2str(i*dt)])
end
%legend('day 0', 'day 1', 'day 2', 'day 3', 'day 4', 'day 5', 'day 6')
%legend boxoff
xlim([ 0 14])
xlabel('infection time (days)')
ylabel('Weighted infectiousness (\beta(\tau)*I(t))')
title('Weighted infectiousness (\beta(\tau)*I(t)) with no isolation')
set(gca,'FontSize',16,'LineWidth',1.5)



%% Test effect of alpha (% removed) on proportion infected over time
trem = 3;
tau_params = horzcat(a,b);
alphavec = [0.8, 0];
alphavec = [0.5:.1:1];
paramsi = params;
for i = 1:length(alphavec)
    paramsi(3) = alphavec(i);
    [y, B, new_inf, R_t, new_exp] = fwd_SEIRD_model(paramsi,tau_params, tvec, y0, dt);
    tot_inf = sum(y(:,3));
    S_t(:,i) = y(:,1);
    ihalfi(i) = find(y(:,1)<0.5*N, 1, 'first');
    thalfi(i) = tvec(ihalf);
end


all_S_t = vertcat(S_t(:,1), S_t(:,2));
lab = vertcat(zeros(length(S_t),1), ones(length(S_t),1));

[b1, logl, H, stats] = coxphfit(lab, all_S_t);
hr = exp(b1)
hrl = exp(b1-stats.se)
hrh = exp(b1+stats.se)

figure;
for i = 1:length(alphavec)
plot(tvec, S_t(:,i)/N, '-', 'LineWidth', 2)
hold on
end
%legend('immediate isolation of symptomatic', 'no isolation', 'Location', 'NorthWest')
legend([num2str(100*alphavec(1)),'% infections removed'], [num2str(100*alphavec(2)),'% infections removed'], [num2str(100*alphavec(3)),'% infections removed'],[num2str(100*alphavec(4)),'% infections removed'],[num2str(100*alphavec(5)),'% infections removed'],[num2str(100*alphavec(6)),'% infections removed'],'Location', 'NorthWest')
legend boxoff
ylabel('Percent not infected')
xlabel('epidemic time (days)')
set(gca,'FontSize',16,'LineWidth',1.5)
%title(['HR for infection = ', num2str(hr),'[' num2str(hrl),',', num2str(hrh),']'])
title('Effect of isolation efficiency')
%% Test effect of more realistic (longer) time to remove

tremvec = [2.5 3.5 4.5 5.5 14];
params = [R_0, gamma, alpha, infend, trem];
tau_params = horzcat(a,b);
paramsi = params;
for i = 1:length(tremvec)
    paramsi(5) = tremvec(i);
    [y, B, new_inf, R_t, new_exp] = fwd_SEIRD_model(paramsi,tau_params, tvec, y0, dt);
    tot_inf = sum(y(:,3));
    S_t(:,i) = y(:,1);
end


all_S_t = vertcat(S_t(:,1), S_t(:,2));
lab = vertcat(zeros(length(S_t),1), ones(length(S_t),1));

[b1, logl, H, stats] = coxphfit(lab, all_S_t);
hr = exp(b1)
hrl = exp(b1-stats.se)
hrh = exp(b1+stats.se)

figure;
for i = 1:length(tremvec)
plot(tvec, S_t(:,i)/N, '-', 'LineWidth', 2)
hold on
end
%legend('immediate isolation of symptomatic', 'no isolation', 'Location', 'NorthWest')
legend('immediate isolation', '1 day to isolate', '2 days to isolate', '3 days to isolate', 'no isolation')
legend boxoff
ylabel('Percent not infected')
xlabel('epidemic time (days)')
set(gca,'FontSize',16,'LineWidth',1.5)
title('Effect of delay in isolation on infection dynamics')
%title(['HR for infection = ', num2str(hr),'[' num2str(hrl),',', num2str(hrh),']'])



%% Ideal vs. Realistic scenario
tremvec = [3, 5];
alphavec = [0.8, 0.5];
paramsi = params;
for i = 1:length(alphavec)
    paramsi(3) = alphavec(i);
    paramsi(5) = tremvec(i);
    [y, B, new_inf, R_t, new_exp] = fwd_SEIRD_model(paramsi,tau_params, tvec, y0, dt);
    tot_inf = sum(y(:,3));
    S_t(:,i) = y(:,1);
    thalf(i)
    pct_infected(i) = 1-(S_t(end,i)./N)
end


all_S_t = vertcat(S_t(:,1), S_t(:,2));
lab = vertcat(zeros(length(S_t),1), ones(length(S_t),1));

[b1, logl, H, stats] = coxphfit(lab, all_S_t);
hr = exp(b1)
hrl = exp(b1-stats.se)
hrh = exp(b1+stats.se)

figure;
for i = 1:length(alphavec)
plot(tvec, S_t(:,i)/N, '-', 'LineWidth', 2)
hold on
end
legend('ideal', 'realistic', 'Location', 'NorthWest')
legend boxoff
ylabel('Percent not infected')
xlabel('epidemic time (days)')
set(gca,'FontSize',16,'LineWidth',1.5)
title(['HR for infection = ', num2str(hr),'[' num2str(hrl),',', num2str(hrh),']'])
title(['Ideal ', num2str(round(100*pct_infected(1),0)),'% infected vs realistic ', num2str(round(100*pct_infected(2),0)),'% infected scenarios'])

%%
figure;
plot(tauvec, w_tau, 'k--', 'LineWidth', 2)
hold on
for i = 2:4:18%:length(tvec)
plot(tauvec, B(i,:), '-', 'LineWidth',2)
text(tauvec(1), B(i,1), ['t=',num2str(i/dt), ' days'])
hold on
end
legend('infectiousness')
xlabel('time since infection (days)')
ylabel('Number of individuals in each stage of infection')
title('Distribution of time since infection')
set(gca,'FontSize',16,'LineWidth',1.5)

figure;
for i = 1:length(tauvec)
plot(tvec, B(:,i), '-')
hold on
end
xlabel('epidemic time (days)')
ylabel('Number of individuals infected')
set(gca,'FontSize',16,'LineWidth',1.5)

figure;
imagesc(B)
colorbar
ylabel('epidemic time (t)(days)')
xlabel('infection time (\tau) (days)')
set(gca,'FontSize',16,'LineWidth',1.5)
