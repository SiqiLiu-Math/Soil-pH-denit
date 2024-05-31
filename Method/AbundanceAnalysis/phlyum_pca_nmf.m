% Do PCA or NMF for phylum level growth and death

clear
load('results/abundance_coarse_grain.mat');

lv = 3;
nn = 1:40;
a00 = Abd_data(lv).Native(:,nn);
a2 = Abd_data(lv).None(:,nn);
a1 = Abd_data(lv).CHL(:,nn);
a0 = zeros(130,length(nn));
for ii=1:10
    kk = (1:13)+(ii-1)*13;
    a0(kk,:) = ones(13,1)*a00(ii,:);
end

lg0 = log(sum(a00,1)+1e-3)/log(10);
wp = transpose(sum(a00,1)+1e-8/10);
lg1 = log(a1+1e-3)/log(10)-log(a0+1e-3)/log(10);
lg2 = log(a2+1e-3)/log(10)-log(a1+1e-3)/log(10);

alg = {'PCA';'PCA';'NMF';'NMF'};
dtp = {'death';'growth';'death';'growth'};
PD = struct('Algorithm',alg,'DataType',dtp,'EigenVector',[],'Weights',[],'Explained',[]);

%% pca
lg1s = lg1([1:65,68:130],:);
lg2s = lg2([1:65,68:130],:);
[H1,W1,~,~,D1,~] = pca(lg1s','Centered',false,'Weights',wp);
[H2,W2,~,~,D2,~] = pca(lg2s','Centered',false,'Weights',wp);
PD(1).Weights = W1';
PD(2).Weights = W2';
PD(1).EigenVector = [H1(1:65,:);NaN(2,40);H1(66:128,:)];
PD(2).EigenVector = [H2(1:65,:);NaN(2,40);H2(66:128,:)];
PD(1).Explained = D1;
PD(2).Explained = D2;

%% nmf
lg01 = max(-lg1,0);
lg02 = max(lg2,0);
[H01,W01,D01] = nnmf(lg01,2,'algorithm','mult','replicates',50);
[H02,W02,D02] = nnmf(lg02,2,'algorithm','mult','replicates',50);
H01(66:67,:) = NaN; H01 = -H01;
H02(66:67,:) = NaN;
PD(3).Weights = W01;
PD(4).Weights = W02;
PD(3).EigenVector = H01;
PD(4).EigenVector = H02;
PD(3).Explained = (1-D01)*100;
PD(4).Explained = (1-D02)*100;

save('results/phylum_pca_nmf_data.mat','lg1','lg2','lg0','ph0','ph1','PD');