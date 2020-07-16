function [y, B, new_infi, expbetai, inf_distrib, Reff] = fwd_SEIRD_model(params,tau_params,tvec, y0, dt, disease)
% This function runs the forward finite difference implementation of the
% SEIRD model, incorporating the time-dependent infectiousness by
% convolving beta(t) with I(t) at each time point...
P = num2cell(params);
Y = num2cell(y0);
[beta0, gamma, alpha, infend, trem]=deal(P{:});
switch disease
    case 'COVID'
    a= tau_params(1);
    b = tau_params(2);
    tind = tvec(1):dt:tvec(end);
    tind = tind';
    w_tau = gampdf(tind, a,b);
    
    case 'SARS'
    mu = tau_params(1);
    sigma = tau_params(2);
    tind = tvec(1):dt:tvec(end);
    tind = tind';
    w_tau = normpdf(tind, mu, sigma);
    
    case 'flu'
    a= tau_params(1);
    b = tau_params(2);
    tind = tvec(1):dt:tvec(end);
    tind = tind';
    w_tau = gampdf(tind, a,b);
    
end
    


beta_tau = beta0*w_tau;

%
[S0, E0, I0, Isr0, R0]=deal(Y{:});


S(1) = S0;
E(1) = E0;
I(1) = I0;
Isr(1) = Isr0;
R(1) = R0;
N = S0+E0+I0+Isr0+R0;




gamma = gamma; % Rate of going from exposed to infectious
beta0 = beta0; % Rate infected spread to susceptible 
delta = 1/(infend-trem); % Rate of recovery from isolated = duration of infection- duration infected not isolated
infend = infend/dt+1; % time to be considered "recovered"
alpha = alpha; % Percent isolated at tremoved
irem = round(trem/dt+1,0); 
% B keeps track of the individuals in each stage of infection (density of 
B=zeros(length(tind), length(tind)); % start with no one infected
B(1,1) = I0; % first time point have one infected at age 1 of infection
%
for t = 2:length(tind)
 

    % Use previus time step to update each compartment
  
        % those in B matrix after time of being infectious become recovered
        recovered = B(t-1, infend);
        % alpha % of those after irem get put into isolation
        isolated = alpha*(B(t-1,irem));
        recoveredisolated = delta*Isr(t-1);
        % remove the isolated from B matrix
        B(t-1,irem) = B(t-1,irem)-(alpha.*B(t-1,irem));
        
         Isr(t) = Isr(t-1) + dt*(isolated-recoveredisolated);
         R(t) = R(t-1) + dt*(recovered + recoveredisolated);
    % Take density of infectious infidivuals over infection age and
    % multiply by weight of infection age
    expbeta = sum(B(t-1,1:end-1)'.*beta_tau(2:end).*dt); % expectation of beta
    pdf_inf =  B(t-1,1:end-1)./sum(B(t-1,1:end-1));
    beta_t = sum(pdf_inf'.*beta_tau(2:end)); % average/expected value of beta 
    % update S, E, and I accordingly 
    new_exposures = expbeta*S(t-1)./(N);
    %new_exposures = expbeta*S(t-1)./(N-Isr(t-1)); % run if you want to remove isoalted from system
    new_inf = gamma*E(t-1);
    S(t) = S(t-1) - dt*(new_exposures);
    E(t) = E(t-1) + dt*(new_exposures - new_inf);
    I(t) = I(t-1) + dt*(new_inf - recovered -isolated);
    % update B and count of new infections
    B(t,1) = new_inf;
    new_infi(t)=new_inf;
    new_exp(t) = new_exposures;
    expbetai(t) = beta_t;
    inf_distrib(:,t) = beta0.*(B(t-1,:).*(w_tau./dt)'); % each column is a distribution of infectionusess
  
    %Re(t) = new_infi(t)/I(t-1);
    % update B by moving each infection along one step form previous
    % infection
    for i = 2:length(tind)
        B(t,i) = B(t-1,i-1); % each time step, push infections along 
    end
    
    
end


y = horzcat(S',E',I', Isr', R');

S(1) = S0;
E(1) = E0;
I(1) = I0;
Isr(1) = Isr0;
R(1) = R0;
N = S0+E0+I0+Isr0+R0;
B2=zeros(length(tind), length(tind)); % start with no one infected
B2(1,1) = I0; % first time point have one infected at age 1 of infection
%
% Estimate Reff
sec_inf=0;
new_inf2(1) = 0;
expbetat(1) = 0;
for t = 2:length(tind)
 

    % Use previus time step to update each compartment
  
        % those in B matrix after time of being infectious become recovered
        recovered = B2(t-1, infend);
        % alpha % of those after irem get put into isolation
        isolated = alpha*(B2(t-1,irem));
        recoveredisolated = delta*Isr(t-1);
        % remove the isolated from B matrix
        B2(t-1,irem) = B2(t-1,irem)-(alpha.*B2(t-1,irem));
        
         Isr(t) = Isr(t-1) + dt*(isolated-recoveredisolated);
         R(t) = R(t-1) + dt*(recovered + recoveredisolated);
    % Take density of infectious infidivuals over infection age and
    % multiply by weight of infection age
    expbeta = sum(B2(t-1,1:end-1)'.*beta_tau(2:end).*dt); % expectation of beta 
    % update S, E, and I accordingly 
    new_exposures = expbeta*S(t-1)./(N);
    expbetat(t) = expbeta;
    %new_exposures = expbeta*S(t-1)./(N-Isr(t-1)); % run if you want to remove isoalted from system
    new_inf = gamma*E(t-1);
    S(t) = S(t-1) - dt*(new_exposures);
    E(t) = E(t-1) + dt*(new_exposures- new_inf);
    I(t) = I(t-1) + dt*(new_inf- recovered - isolated);
    % update B and count of new infections
    new_inf2(t) = new_inf;
    sec_inf = sec_inf + new_inf;
    %Re(t) = new_infi(t)/I(t-1);
    % update B by moving each infection along one step form previous
    % infection
    for i = 2:length(tind)
        B2(t,i) = B2(t-1,i-1); % each time step, push infections along 
    end
    B2(t,1) = 0;
    
end

Reff = sec_inf./I0;


end
