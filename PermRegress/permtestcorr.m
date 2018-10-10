function [pval,grotperms,side] = permtestcorr(Xin,Yin,Nperm,confounds,cs,mode)
% tests if the correlation between Xin and Yin is different from zero
% mode == 1 to base on sign flipping
% mode == 2 to base on permutations
% Diego Vidaurre, University of Oxford (2015)

N = length(Yin);
if (nargin>3) && ~isempty(confounds)
    confounds = confounds - repmat(mean(confounds),N,1);
    %Xin = Xin - confounds * pinv(confounds) * Xin;
    Yin = Yin - confounds * pinv(confounds) * Yin;
end
if (nargin<5)
    cs = [];
else
    if ~isempty(cs)
        if size(cs,2)>1 % matrix format
            [allcs(:,2),allcs(:,1)]=ind2sub([length(cs) length(cs)],find(cs>0));
            [grotMZi(:,2),grotMZi(:,1)]=ind2sub([length(cs) length(cs)],find(tril(cs,1)==1));
            [grotDZi(:,2),grotDZi(:,1)]=ind2sub([length(cs) length(cs)],find(tril(cs,1)==2));
        else
            allcs = [];
            nz = cs>0;
            gr = unique(cs(nz)); ngr = length(gr);
            for g=gr'
                ss = find(cs==g);
                for s1=ss
                    for s2=ss
                        allcs = [allcs; [s1 s2]];
                    end
                end
            end
        end
    end
end
if (nargin<6)
    mode = 2;
end

Xin = Xin - mean(Xin); Xin = Xin / std(Xin);
Yin = Yin - mean(Yin); Yin = Yin / std(Yin);

if mode == 1
    [pval,~,grotperms] = permtestmean(Yin.*Xin,Nperm);
else
    grotperms = zeros(Nperm,1);
    for perm=1:Nperm
        Yin0 = Yin;
        if (perm>1)
            if isempty(cs)           % simple full permutation with no correlation structure
                rperm = randperm(N);
                Yin=Yin0(rperm,:);
            else
                PERM=zeros(1,N);
                perm1=randperm(size(grotMZi,1));
                for ipe=1:length(perm1) % you basically permute within families
                    if rand<0.5, wt=[1 2]; else wt=[2 1]; end
                    PERM(grotMZi(ipe,1))=grotMZi(perm1(ipe),wt(1));
                    PERM(grotMZi(ipe,2))=grotMZi(perm1(ipe),wt(2));
                end
                perm1=randperm(size(grotDZi,1));
                for ipe=1:length(perm1)
                    if rand<0.5, wt=[1 2]; else wt=[2 1]; end
                    PERM(grotDZi(ipe,1))=grotDZi(perm1(ipe),wt(1));
                    PERM(grotDZi(ipe,2))=grotDZi(perm1(ipe),wt(2));
                end
                from=find(PERM==0);  pto=randperm(length(from));  to=from(pto);  PERM(from)=to;
                Yin=Yin0(PERM,:);
            end
        end
        grotperms(perm) = abs(corr(Xin,Yin));
        if perm==1, side = sign(corr(Xin,Yin)); end
    end
    pval = sum(grotperms>=grotperms(1)) / (Nperm+1);
end

end