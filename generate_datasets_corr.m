function [Y,x,beta] = generate_datasets_corr(N,T,R,P,maxbeta,sigma,eta,nu)
% The general setting is:
%
% - Hidden variable beta, (N by 1), containing a certain modulating parameter 
%       (say a certain neural variable of interest, e.g. FC)
% - Observed value x (N by 1), which is a variable possibly related to beta
%       (say this is a behavioural variable which presumably relates to FC, e.g. IQ)
% - Hidden variable S_i, (T_i by Q), corresponding to certain neural process 
%       of interest, for each subject i=1...N, which is modulated by beta_i
%       through some unknown function f()
%       (for example S_i is the brain process under beta_i's FC) 
% - Observed data, D_i (R_i by Q), a noisy version of S_i (say imaging data)
% - Estimated values Y, (N by P), is derived as Y_ij = g(D_i), i=1,...,N, 
%       which we can repeat a number of times j=1...P
%       (say a measure of FC estimated from noisy brain imaging data)
%
% We don't have access to neither S, beta or f, but we have a mechanism to
% estimate Y_i = g(D). Function g() is cheap to evaluate, but is noisy
% and returns a different estimation D each time, i.e. 
% it has a certain variance, Var(g(S)) = epsilon.
% 
% The goal is then to test for the relationship between y and beta,
% leaning on the dependence between (observable) Y and beta.
% Alternative hypothesis: x (IQ) and beta (FC) are related
% In the example here, this is true, but we can fail to reject the null 
% if the chain of dependences Y<-D<-S<-beta is weak 
% 
% This function generates data following this setting:
% - beta is a scaling coefficient, uniformly sampled from [-maxbeta,maxbeta]
% - x is the (noisy corr coeff) computed as y = h(beta) + Gaussian_noise(0,sigma),
%       such that h(beta) = beta / sqrt(beta^2 + eta^2), 
%       which is the analytical function that maps from a
%       scaling coefficient to a correlation coefficient
% - D_i is (T by 2) is latent data 
%   . D_i(:,1) is white noise (mean 0, std deviation 1)
%   . D_i(:,2) = D_i(:,1) * beta(i) + eta * white noise
% - Dhat_i is a noisy subsample (R samples) of D,
%       such that Dhat_i = subsample(D_i,R) + Gaussian_noise(0,nu)
% - Y is (N by 1), such that Y(i) = corr(Dhat_i(:,1),Dhat_i(:,2))
%
% Then, in summary
% - eta and maxbeta control the correlation coefficient 
% - sigma controls how noisy is the observation of the correlation coefficient
% - nu and R control how noisy is our sampling mechanism 
% - P is how many noisy measures we want to have
%
% (Note here we have a noisy subsample of the data each time; in the
% typical HMM type-of-analysis scenario it is the estimation of Y what is
% noisy and done various times. But there is no loss of generality as the
% noisy estimation happens in the Y<-Dhat<-D<-beta chain anyway)
%
% Diego Vidaurre (2018)

Y = zeros(N,P);
beta = rand(N,1) * 2*maxbeta - maxbeta; 
x = beta ./ sqrt(beta.^2 + eta.^2) + sigma * randn(N,1);
% create subsampling
r = zeros(R,P);
for j = 1:P, r(:,j) = randsample(N,R); end
for i = 1:N
   S = zeros(T,2); 
   S(:,1) = randn(T,1);
   S(:,2) = S(:,1) * beta(i) + eta * randn(T,1); 
   for j = 1:P
       D = S(r(:,j),:) + nu * randn(R,2);
       Y(i,j) = atanh(corr(D(:,1),D(:,2)));
   end
end

end






