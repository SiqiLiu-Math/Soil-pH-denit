% Load nitrate dynamics data from matlab table to matlab cell

clear
load('../../Data/ProcessedData/Denitrification_data_20soil.mat');
sample_size = size(DN_none,1);

fdata = cell(sample_size,2); % CHL-/+ data

for nn=1:sample_size
    t1 = cell2mat(table2array(DN_none(nn,'time_none')));
    a1 = cell2mat(table2array(DN_none(nn,'no3_none')));
    t2 = cell2mat(table2array(DN_chl(nn,'time_chl')));
    a2 = cell2mat(table2array(DN_chl(nn,'no3_chl')));
    fdata{nn,1} = [t1;a1];
    fdata{nn,2} = [t2;a2];
end

save('results/nitrate_data_for_fit.mat','fdata');