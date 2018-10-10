function pvals = permtestmass_regress(Yin,Xin,Nperm,Perms,conf,verbose)
% 
% It tests all replications (column) of Xin vs each variable (column) of Yin, 
% using linear regression embedded in permutation testing 
%
% ARGUMENTS
% - Yin: Noisy replications
% - Xin: Observed behavioural variables.
% - Nperm is the number of permutations
% - Perms: precomputed permutations, (Nsubjects by Nperm)
% - conf are the confounds
%
% - pvals contains 1 p-value per column in Xin;
%
% Diego Vidaurre, University of Oxford (2017)


if nargin<4 
    Perms = [];
end
if nargin<5 
    conf = [];
end
if nargin<6
    verbose = 0;
end

N = size(Xin,1);
R = size(Yin,2);
P = size(Xin,2);

if ~isempty(conf)
    conf = conf - repmat(mean(conf),N,1);
    Yin = Yin - repmat(mean(Yin),N,1);
    Yin = Yin - conf * pinv(conf) * Yin;
    %if any(~is_binary)
    %    Yin(:,~is_binary) = Yin(:,~is_binary) - repmat(mean(Yin(:,~is_binary)),N,1);
    %    Yin(:,~is_binary) = Yin(:,~is_binary) - conf * pinv(conf) * Yin(:,~is_binary);
    %end
end

Yin = zscore(Yin);
Xin = zscore(Xin);
sqerr = zeros(Nperm,P);
% We should do ridge regression instead
%if R > N, [~,Yin] = pca(Yin,'NumComponents',round(N/2)); end 
alpha = 0;
if R > N, alpha = 0.001; end

gram = (Yin' * Yin + alpha * eye(R)) \ Yin';

% Unpermuted
beta = gram * Xin;
sqerr(1,:) = sum((Xin - Yin * beta).^2);
if verbose, disp('1'); end

% Outer perm loop, rest of permutations
X = Xin; 
parfor i = 2:Nperm
    if isempty(Perms)
        perm = randperm(N);
    else
        perm = Perms(:,randperm(size(Perms,1),1));
    end
    Xin = X(perm,:);
    beta = gram * Xin;
    sqerr(i,:) = sum((Xin - Yin * beta).^2);
    if verbose, disp(num2str(i1)); end
end

pvals = sum( repmat(sqerr(1,:),Nperm,1) >= sqerr ) ./ (1 + Nperm) ;

end


