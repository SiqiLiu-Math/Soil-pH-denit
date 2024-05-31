% compute error distribution of fitted results

clear
load('results/model_parameters_fitting_nitrate.mat');
load('results/nitrate_data_for_fit_cleaned.mat');
load('../../Data/ProcessedData/Denitrification_data_20soil.mat');
sample_size = length(paras);

ers = cell(sample_size,1);
for ii=1:sample_size
    er = er_function(paras(ii,:),fdata{ii,1},fdata{ii,2});
    ers{ii} = er';
end
error_per_datapoint = cell2mat(ers);

error_per_condition = zeros(sample_size/3,1);
for ii=1:sample_size/3
    kk = ii*3-[2,1,0];
    erc = cell2mat(ers(kk,1));
    error_per_condition(ii) = mean(erc.^2);
end
ph_per_condition = table2array(DN_none(1:3:end,3:4));

save('results/error_distribution.mat','error_per_datapoint','error_per_condition','ph_per_condition');

%%%%%%%%%%%
function er = er_function(paras,fd1,fd2)

a0c = paras(1);
a0n = paras(2);
x0 = paras(3);
ts = paras(4);
ga = 4.8;
t1 = fd1(1,:);
a1 = fd1(2,:);
t2 = fd2(1,:);
a2 = fd2(2,:);

A1 = a0n-x0/ga.*(exp(ga*min(t1,ts))-1)-x0*exp(ga*ts)*max(t1-ts,0);
A1 = max(A1,0);
A2 = max(a0c-x0*t2,0);

er = [A1-a1,A2-a2];

end