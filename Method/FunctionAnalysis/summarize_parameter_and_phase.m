% Summarize the model parameters and name phases for each condition

clear
load('results/model_parameters_fitting_nitrate.mat')
load('../../Data/ProcessedData/Denitrification_data_20soil.mat');

ll = length(paras)/3;
x0 = reshape(table2array(model_parameters(:,3)),[3,ll]);
c0 = reshape(table2array(model_parameters(:,4)),[3,ll]);
x0 = transpose(median(x0,1));
c0 = transpose(median(c0,1));

phase = 2*ones(ll,1);
phase(x0<0.05) = 1;
phase(c0>0.3) = 3;
lgx0 = log(x0+2.4e-2)/log(10);
lgc0 = log(c0+1e-2)/log(10);

ph0 = table2array(DN_none(1:3:end,'ph_soil'));
ph1 = table2array(DN_none(1:3:end,'ph_none'));
unt = table2array(DN_none(1:3:end,'unit'));
sid = table2array(DN_none(1:3:end,'soil_id'));

save('results/function_parameter_and_phases.mat','sid','ph0','ph1','unt','x0','c0','lgx0','lgc0','phase');