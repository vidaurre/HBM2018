% Script to generate the data for Figure 1. 

R = 1000; 
T = 100; 
c = [0 0.1 0.2]; L = length(c);
tests = zeros(R,L);
tests_fdr = zeros(R,L);
test_nested = zeros(1,L);
empirical_corr = zeros(R,L);

addpath('PermRegress/')


for i = 1:L
    X = randn(T,R);
    y = randn(T,1);
    ci = c(i);
    for j = 1:R
        if ci>0, X(:,j) = X(:,j) + ci * rand * y; end
        empirical_corr(j,i) = corr(X(:,j),y);
        tests(j,i) = permtestcorr(X(:,j),y,10000);
        if rem(j,50)==0, disp([num2str(i) ' ' num2str(j)]); end
    end
    tests_fdr(:,i) = mafdr(tests(:,i),'BHFDR',1);
    test_nested(i) = permtestmass_NPC(X,y,10000);
    i
end
    
save('out/demo.mat','empirical_corr','tests','tests_fdr','test_nested')
%% 

load('out/demo.mat')

figure(1)
for i = 1:L
    subplot(2,L,i); [b,x]=hist(empirical_corr(:,i)); 
    bar(x,b,'FaceColor','r')
    set(gca,'FontSize',17); ylim([0 1.1*max(b)]); %'YTIck',[]
    subplot(2,L,i+L); [b,x]=hist(tests(:,i)); 
    bar(x,b,'FaceColor','r')
    set(gca,'FontSize',17); ylim([0 1.1*max(b)]); xlim([0 1]); %YTIck',[]
end

%%

mean(tests)
mean(tests_fdr)

exp(mean(log(tests)))
exp(mean(log(tests_fdr)))


