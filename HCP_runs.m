replications = 100;

DirData = '~/MATLAB/data/HCP/TimeSeries/3T_HCP820_MSMSulc_d50_ts2/';
DirHCP = '~/MATLAB/data/HCP/';
DirOut = '/home/diegov/MATLAB/HCP/HMM/out/multiple_repetitions/';

addpath(genpath('~/MATLAB/HMM-MAR'))

% loading data
files = dir(DirData);
files = files(3:end); 
N = length(files);
f = cell(N,1);
T = cell(N,1);
for ifile = 1:N
    f{ifile} = strcat(DirData,files(ifile).name);
    T{ifile} = [1200 1200 1200 1200];
end 

% HMM options - for specific info check
% https://github.com/OHBA-analysis/HMM-MAR/wiki/User-Guide 
options = struct();
% Specific Stochastic HMM options
options.BIGmincyc = 10;
options.BIGcyc = 100;
options.BIGundertol_tostop = 1;
options.BIGNbatch = round(N/30);
options.BIGtol = 1e-7;
options.BIGcyc = 200;
options.BIGundertol_tostop = 5;
options.BIGdelay = 1;
options.BIGforgetrate = 0.7;
options.BIGbase_weights = 0.9;
options.BIGverbose = 0;
% Training HMM options
options.K = 12; 
options.covtype = 'full';
options.order = 0;
options.zeromean = 0;
options.decodeGamma = 0;
options.tol = 1e-7;
options.cyc = 10;
options.DirichletDiag = 10;
options.initcyc = 0; % this indicates a random init, which we normally wouldn't do
options.initrep = 0; % this indicates a random init, which we normally wouldn't do
options.verbose = 0; 
 
% Different subdirs for different parametrisations
if options.zeromean==1 && strcmp(options.covtype,'full')
    subdir = [subdir 'zeromean1_full/'];
elseif options.zeromean==0 && strcmp(options.covtype,'uniquefull')
    subdir = [subdir 'zeromean0_uniquefull/'];
else
    subdir = [subdir 'zeromean0_full/'];
end

% for monitoring purposes, check free energy
FreeEnergy = zeros(replications,1);

% fractional occupancy, and switching rate
FO = zeros(820,12);
dFO = zeros(820,12);

% Run the HMM a lot of times
for r = 1:replications
    [hmm,Gamma] = hmmmar(f,T,options);
    FO = zeros(820,12);
    for j = 1:820
        ind = (1:4800) + (j-1) * 4800;
        FO(j,:) = mean(Gamma(ind,:));
        dFO(j,:) = mean(abs(Gamma(ind(2:end),:) - Gamma(ind(1:end-1),:)));
    end
    save([DirOut subdir 'run_' num2str(r) '.mat'],'FO','dFO','hmm')
    disp(num2str(r))
end