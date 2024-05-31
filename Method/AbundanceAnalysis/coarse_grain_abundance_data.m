% Coarse-grain the abundance data according to taxa (phylum to species)

clear
load('results/abundance_sort_by_native_taxa.mat');
noc = size(taxa_id,1);

txlevel = {'Biomass';'Kingdom';'Phylum';'Class';'Order';'Family';'Genus';'Species'};
Abd_data = struct('Level',txlevel,'Name',[],'Native',[],'None',[],'CHL',[],'Tree',[]);
tid = [ones(noc,1),taxa_id];

for ii=1:8
    [~,uu] = unique(tid(:,ii));
    Abd_data(ii).Name = table2cell(taxa_inf(uu,2+ii));
    abd_native = zeros(10,length(uu));
    abd_none = zeros(130,length(uu));
    abd_chl = zeros(130,length(uu));
    for jj=1:noc
        kk = tid(jj,ii);
        abd_native(:,kk) = abd_native(:,kk)+abundance_native(:,jj);
        abd_none(:,kk) = abd_none(:,kk)+abundance_none(:,jj);
        abd_chl(:,kk) = abd_chl(:,kk)+abundance_chl(:,jj);
    end
    Abd_data(ii).Native = abd_native;
    Abd_data(ii).None = abd_none;
    Abd_data(ii).CHL = abd_chl;
    Abd_data(ii).Tree = tid(uu,1:ii);
end
Abd_data(1).Name = 'Biomass';

save('results/abundance_coarse_grain.mat','Abd_data','ph1','ph0');
