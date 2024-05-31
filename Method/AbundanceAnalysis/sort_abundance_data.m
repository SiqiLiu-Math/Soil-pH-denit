% Decompose abundance matrix into: 1) none; 2) chl; 3) native.
% Combine three replicates. Rescale by native. Sort taxa by native.

clear
load('../../Data/ProcessedData/Abundance_data_10soil.mat');


num_of_column = size(abundance_matrix,2);
column_id = table2array(column_info(:,'Column_Id'));
abundance_none = zeros(130,num_of_column);
abundance_chl = zeros(130,num_of_column);
abundance_native = zeros(10,num_of_column);
ph0 = zeros(130,1);
ph1 = zeros(130,1);
for ii=1:10
    k0 = (85:87)+(ii-1)*87;
    k1 = (1:39)+(ii-1)*87;
    k2 = (40:78)+(ii-1)*87;
    sc = sum(abundance_matrix(k0,:),'all');
    abundance_native(ii,:) = sum(abundance_matrix(k0,:),1)/sc;
    for jj=1:13
        kk = (ii-1)*13+jj;
        l1 = k1(jj*3-[2,1,0]);
        l2 = k2(jj*3-[2,1,0]);
        ph0(kk) = table2array(row_info(l1(1),'Soil_pH'));
        ph1(kk) = table2array(row_info(l1(1),'Perturbed_pH'));
        abundance_none(kk,:) = sum(abundance_matrix(l1,:),1)/sc;
        abundance_chl(kk,:) = sum(abundance_matrix(l2,:),1)/sc;
        if kk==100
            l2 = [l2(1),l2(3)];
            abundance_chl(kk,:) = sum(abundance_matrix(l2,:),1)/sc*3/2;
        end
    end
end


txid = table2array(column_info(:,12:18));
asv_abundance = transpose(sum(abundance_native));
sum_abd = cell(7,1);
for jj=1:7
    mm = txid(end,jj);
    sum_abd{jj} = zeros(mm,1);
end
for ii=1:num_of_column
    for jj=1:7
        kk = txid(ii,jj);
        sum_abd{jj}(kk) = sum_abd{jj}(kk)+asv_abundance(ii);
    end
end
Species_abundance = sum_abd{7}(txid(:,7));
Genus_abundance = sum_abd{6}(txid(:,6));
Family_abundance = sum_abd{5}(txid(:,5));
Order_abundance = sum_abd{4}(txid(:,4));
Class_abundance = sum_abd{3}(txid(:,3));
Phylum_abundance = sum_abd{2}(txid(:,2));
Kingdom_abundance = sum_abd{1}(txid(:,1));
taxa_abd = table(column_id,Kingdom_abundance,Phylum_abundance,Class_abundance,Order_abundance,Family_abundance,Genus_abundance,Species_abundance,asv_abundance);
taxa_abd = sortrows(taxa_abd,"asv_abundance",'descend');
taxa_abd = sortrows(taxa_abd,"Species_abundance",'descend');
taxa_abd = sortrows(taxa_abd,"Genus_abundance",'descend');
taxa_abd = sortrows(taxa_abd,"Family_abundance",'descend');
taxa_abd = sortrows(taxa_abd,"Order_abundance",'descend');
taxa_abd = sortrows(taxa_abd,"Class_abundance",'descend');
taxa_abd = sortrows(taxa_abd,"Phylum_abundance",'descend');
taxa_abd = sortrows(taxa_abd,"Kingdom_abundance",'descend');
new_id = table2array(taxa_abd(:,1));
taxa_inf = column_info(new_id,1:11);
abundance_none = abundance_none(:,new_id);
abundance_chl = abundance_chl(:,new_id);
abundance_native = abundance_native(:,new_id);
txid = txid(new_id,:);
taxa_id = ones(size(txid));
for ii=2:size(txid,1)
    df = (txid(ii,:)~=txid(ii-1,:));
    taxa_id(ii,:) = taxa_id(ii-1,:)+df;
end

save('results/abundance_sort_by_native_taxa.mat','taxa_id','taxa_inf','abundance_native','abundance_chl','abundance_none','taxa_abd','ph0','ph1');
