addpath('PermRegress')

Ns = [50 200 1000]; % No. of subjects
T = 10000; % Original time
R = 20; % Time of Dhat
P = 100; % number of replications
sigma = 0.5; % noise in the relation between beta and y
eta =  1.0; % parameter modulating FC (eta and maxbeta define correlation)
%           y = beta / (beta^2 + eta^2)
% nu, with R, controls the noisiness of the replicating 
%           (the lower, the better we do)

%% Playing with nu


maxbeta = 0.2; % range of the regress coefficients
repetitions = 10;
nu_values = 0.25:0.05:1.5; Nnu = length(nu_values);

for N = Ns
    
    AllPvals = zeros(repetitions,Nnu,P);
    Pvals0 = zeros(repetitions,Nnu);
    Pvals = zeros(repetitions,Nnu);
    PvalsRegression = zeros(repetitions,Nnu);
    corr_x_beta = zeros(repetitions,Nnu);
    corr_Y_x = zeros(repetitions,Nnu,P);
    corr_Y_beta = zeros(repetitions,Nnu,P);
    
    for rep = 1:repetitions
        for inu = 1:Nnu
            nu = nu_values(inu);
            [Y,x,beta] = generate_datasets_corr(N,T,R,P,maxbeta,sigma,eta,nu);
            corr_x_beta(rep,inu) = corr(x,beta);
            corr_Y_x(rep,inu,:) = corr(Y,x);
            corr_Y_beta(rep,inu,:) = corr(Y,beta);
            [pv,~,pv0] = permtestmass_NPC(Y,x,10000);
            AllPvals(rep,inu,:) = pv0;
            Pvals(rep,inu) = pv;
            Pvals0(rep,inu) = mean(pv0);
            pvr = permtest_regress(Y,x,10000);
            PvalsRegression(rep,inu) = pvr; 
            %disp([num2str(rep) ' : ' num2str(inu)])
        end
    end
    
    fname = ['out/synthetic_runs_N' num2str(N) '.mat'];
    
    save(fname,'Pvals','Pvals0',...
       'corr_x_beta','corr_Y_x','corr_Y_beta','AllPvals','PvalsRegression')
    
    %save(fname,'PvalsRegression','-append')
    
    disp(['Done with N=' num2str(N) ])
    
end

%% Same simulation but with absolute no relationship

maxbeta = 0; % range of the regress coefficients
repetitions = 10;
nu_values = 0.25:0.05:1.5; Nnu = length(nu_values);

for N = 200
        
    AllPvals = zeros(repetitions,Nnu,P);
    Pvals0 = zeros(repetitions,Nnu);
    Pvals = zeros(repetitions,Nnu);
    PvalsRegression = zeros(repetitions,Nnu);
    corr_y_beta = zeros(repetitions,Nnu);
    corr_X_y = zeros(repetitions,Nnu,P);
    corr_X_beta = zeros(repetitions,Nnu,P);
    
    for rep = 1:repetitions
        for inu = 1:Nnu
            nu = nu_values(inu);
            [Y,x,beta] = generate_datasets_corr(N,T,R,P,maxbeta,sigma,eta,nu);
            corr_x_beta(rep,inu) = corr(x,beta);
            corr_Y_x(rep,inu,:) = corr(Y,x);
            corr_Y_beta(rep,inu,:) = corr(Y,beta);
            [pv,~,pv0] = permtestmass_NPC(Y,x,10000);
            AllPvals(rep,inu,:) = pv0;
            Pvals(rep,inu) = pv;
            Pvals0(rep,inu) = mean(pv0);
            pvr = permtest_regress(Y,x,10000);
            PvalsRegression(rep,inu) = pvr; 
            %disp([num2str(rep) ' : ' num2str(inu)])
        end
    end
    
    fname = ['out/synthetic_runs_0_N' num2str(N) '.mat'];
    
    save(fname,'Pvals','Pvals0',...
        'corr_x_beta','corr_Y_x','corr_Y_beta','AllPvals','PvalsRegression')
    
    %save(fname,'PvalsRegression','-append')
    
end

