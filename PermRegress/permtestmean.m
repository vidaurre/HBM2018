function [pval,my,grotperms] = permtestmean(Yin,Nperm,confounds)
% tests if Yin has mean higher (or lower) than zero. 
% Diego Vidaurre, University of Oxford (2015)

N = length(Yin);
my = mean(Yin);
if (nargin>2) && ~isempty(confounds)
    confounds = confounds - repmat(mean(confounds),N,1);
    Yin = Yin - confounds * pinv(confounds) * Yin;
    my = mean(Yin); 
else
    confounds=[];
end

grotperms = zeros(Nperm,1);
grotperms(1) = my; 

for iperm = 2:Nperm
    Y = Yin; 
    flip = binornd(1,0.5,N,1)==1;
    Y(flip) = - Y(flip);
    grotperms(iperm) = mean(Y);
end

pval = zeros(1,2);
pval(1) = sum(grotperms<=grotperms(1)) / (Nperm+1);
pval(2) = sum(grotperms>=grotperms(1)) / (Nperm+1);
pval = min(pval);

end