HCPdatadir = '../HCP/'; % set your pathways here
addpath('PermRegress')

inject_noise = 0; % if 1, half of the replications will be randomly generated

pvals = cell(4,1);   
pvals_rows = cell(4,1);   
pvals_cols = cell(4,1);   
Pvals_all = cell(4,1);   
pvals_dec = cell(4,1);   
pvals_rows_dec = cell(4,1);   
pvals_cols_dec = cell(4,1);   
Pvals_all_dec = cell(4,1);   

pvals_cols_regr = cell(4,1);   
pvals_cols_dec_regr = cell(4,1);  

if inject_noise, strsub = '_perturbed';
else, strsub = '';
end

% focus on covariance-based HMMs
DirOut = '~/MATLAB/HCP/HMM/out/multiple_repetitions/zeromean0_full/';
nruns = 100;
Y = [];
grouping = [];
D = []; % state switching
for j = 1:nruns
    load([DirOut 'run_' num2str(j) '.mat'])
    K = size(FO,2);
    Y = [Y FO];
    grouping = [grouping j*ones(1,K)];
    D = [D sum(dFO,2)];
end
N = size(Y,1);
if inject_noise > 0 % random runs
    for j = round(nruns/2)+1:nruns
        D(:,j) = D(randperm(size(D,1),size(D,1)),j);
    end
end

% load behavioural variables 
vars = dlmread([HCPdatadir 'scripts900/vars.txt'],' ');
% load twin information  
twins = dlmread([HCPdatadir 'scripts900/twins.txt'],' ');
twins = twins(2:end,2:end);
% headers info
fid = fopen([HCPdatadir 'column_headers.txt']);
headers = textscan(fid,'%s');
headers = headers{1}; headers(1) = []; % because it takes the 1st row in two elements
% set of permutations created with PALM
load([HCPdatadir 'Permutations/Pset_full.mat']);

Drop = false(N,1);
NoFamilyInfo = [108525, 116322, 146331, 168240, 256540,...
    280941, 524135, 604537, 657659, 867468];
for n = NoFamilyInfo, Drop(vars(:,1)==n) = true; end
vars = vars(~Drop,:);
twins = twins(~Drop,~Drop);
Y = Y(~Drop,:); D = D(~Drop,:);
confounds = [3 4 7]; % sex, age, motion
conf = vars(:,confounds); 
index_vars = 3:264;
valid = mean(~isnan(vars(:,index_vars))); % contain no NAN
index_vars = index_vars(valid>0.75);

X = vars(:,index_vars);
keep = sum(isnan(X),2)==0 & ~isnan(sum(conf,2));
keep_v = var(X(keep,:)) > 0;
X = X(:,keep_v);
index_vars = index_vars(keep_v);
headers = headers(index_vars);
permutations = removeSubjfromPerm(Pset_full,~keep);

% deconfound motion
motion = vars(keep,7);
motion = motion - mean(motion);
Xdec = X(keep,index_vars~=7);
Xdec = Xdec - repmat(mean(Xdec),size(Xdec,1),1);
b = (motion' * Xdec) / (motion' * motion);
Xdec = Xdec - motion * b;
headers_dec = headers(index_vars~=7);

% without deconfounding
[pvals_cols,pvals_cols_FWR,Pvals_all] = ...
    permtestmass_NPC(D(keep,:),X(keep,:),10000,permutations);
% in deconfounded space
[pvals_cols_dec,pvals_cols_FWR_dec,Pvals_all_dec] = ...
    permtestmass_NPC(D(keep,:),Xdec,10000,permutations);

% using regression-based permutation testing
pvals_cols_regr = permtest_regress(D(keep,:),X(keep,:),10000,permutations);
pvals_cols_dec_regr = permtest_regress(D(keep,:),Xdec,10000,permutations);


save(['out/run_permtestmassr' strsub '.mat'],...
    'pvals_cols','pvals_cols_FWR','Pvals_all',...
    'pvals_cols_dec','pvals_cols_FWR_dec','Pvals_all_dec',...
    'pvals_cols_regr','pvals_cols_dec_regr',...
    'headers_dec','headers')

