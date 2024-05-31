% Combine function data and phylum data

clear
load('results/abundance_coarse_grain.mat');
load('../Function_fitting_double_cr_model/results/model_parameters_fitting_nitrate.mat');
load('../../Data/ProcessedData/Abundance_data_10soil.mat');


sid = table2array(row_info(1:87:end,'Soil_Id'));
fdat = zeros(130,2);
for ii=1:length(sid)
    ss = sid(ii);
    for jj=1:13
        kk = (ss-1)*39+(jj-1)*3+[1,2,3];
        ll = (ii-1)*13+jj;
        fd = table2array(model_parameters(kk,3:4));
        fdat(ll,:) = median(fd);
    end
end
phase = 2*ones(130,1);
phase(fdat(:,1)<0.05) = 1;
phase(fdat(:,2)>0.3) = 3;
lgx0 = log(fdat(:,1)+2.4e-2)/log(10);
lgc0 = log(fdat(:,2)+1e-2)/log(10);

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
lgn1 = log(a1(:,2)+a1(:,9)+1e-3)/log(10)-log(a0(:,2)+a0(:,9)+1e-3)/log(10);
lgn2 = log(a2(:,2)+a2(:,9)+1e-3)/log(10)-log(a1(:,2)+a1(:,9)+1e-3)/log(10);
lgm1 = log(a1(:,4)+1e-3)/log(10)-log(a0(:,4)+1e-3)/log(10);
lgm2 = log(a2(:,4)+1e-3)/log(10)-log(a1(:,4)+1e-3)/log(10);

save('results/function_and_phylum_data.mat','lgn1','lgn2','lgm1','lgm2','lgc0','lgx0','phase','ph0','ph1','fdat');