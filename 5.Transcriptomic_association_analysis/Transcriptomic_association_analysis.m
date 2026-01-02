
% This code acquires significantly related genes with MFC using PLS 

ROI = {'A1C'
'DFC'
'IPC'
'ITC'
'M1C'
'MFC'
'OFC'
'S1C'
'STC'
'V1C'
'VFC'};

Gene_expr_pre = readmatrix('./data/gene_expression.csv');
load('./data/Gene_label.mat')
load('./data/MFC.mat')

genes=Gene_label; % this needs to be imported first
geneindex=1:length(Gene_label);
bootnum=5000;
X=log2(Gene_expr_pre'+1);
Y=zscore(MFC);
dim=10;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);  %cumsum(PCTVAR,2)
for i =1:length(Gene_label)
    XL_new(i) = corr(X(:,i),XS(:,1));
end

%store regions' IDs and weights in descending order of weight for both components:
[r1(1,1),~]=corr(XS(:,1),Y);
[r1(2,1),~]=corr(XS(:,2),Y);
%align PLS components with desired direction for interpretability 
if r1(1,1)<0
    stats.W(:,1)=-1*stats.W(:,1);
    XS(:,1)=-1*XS(:,1);
end
if r1(2,1)<0
    stats.W(:,2)=-1*stats.W(:,2);
    XS(:,2)=-1*XS(:,2);
end

[PLS1w,x1] = sort(stats.W(:,1),'descend');
PLS1ids=genes(x1);
geneindex1=geneindex(x1);
[PLS2w,x2] = sort(stats.W(:,2),'descend');
PLS2ids=genes(x2);
geneindex2=geneindex(x2);
XL_new = XL_new(x1);

%define variables for storing the (ordered) weights from all bootstrap runs
PLS1weights=[];
PLS2weights=[];

% Perm test - PLS1 PCTVAR
parfor i=1:bootnum
    
    Y_perm = Y(randperm(size(Y,1)), :);
    [~,~,~,~,~,PCTVAR,~,~]=plsregress(X,Y_perm,dim); %perform PLS for resampled data
      
    PCTV(i,1) = PCTVAR(1,1);
       
end

% start bootstrap
parfor i=1:bootnum
    
    myresample = randsample(size(X,1),size(X,1),1);
    res(i,:)=myresample; %store resampling out of interest
    Xr=X(myresample,:); % define X for resampled subjects
    Yr=Y(myresample,:); % define X for resampled subjects
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(Xr,Yr,dim); %perform PLS for resampled data
      
    temp=stats.W(:,1);%extract PLS1 weights
    newW=temp(x1); %order the newly obtained weights the same way as initial PLS 
    if corr(PLS1w,newW)< 0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW=-1*newW;
    end
    PLS1weights(:,i)=newW;%store (ordered) weights from this bootstrap run
    
end

%get standard deviation of weights from bootstrap runs
PLS1sw=std(PLS1weights');

%get bootstrap weights
temp1=PLS1w./PLS1sw'; 
p1 = 2 * (1 - normcdf(abs(temp1)));   % p-values for PLS1 weights
%order bootstrap weights (Z) and names of terms
[Z1 ind1]=sort(temp1,'descend');
P1 = p1(ind1);
PLS1=PLS1ids(ind1); % ordered terms ID
XL_new = XL_new(ind1);  % calculated via correlation between X_orignal and XScore
fdr = mafdr(P1,'BHFDR',1);

index_sig = find(fdr<0.001 & abs(Z1)>3);
orderedtable_PLS = table(PLS1(index_sig),P1(index_sig),fdr(index_sig),Z1(index_sig),XL_new(index_sig)','VariableNames',{'GenesID'; 'P_values'; 'P_fdr';'Loadings_Zscore';'Loadings_corr'});
writetable(orderedtable_PLS(:,[1,4]),'./result/OrderedTerms_PLS_Loadings_Pfdr001_prenatal.csv')


