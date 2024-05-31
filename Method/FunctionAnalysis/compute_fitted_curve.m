% compute the curves of fitted dynamics

clear
load('results/model_parameters_fitting_nitrate.mat')
load('results/nitrate_data_for_fit_cleaned.mat')
num_of_rows = size(paras,1);

T = 0:0.001:4;
a1s = zeros(num_of_rows,length(T));
a2s = zeros(num_of_rows,length(T));
for ii=1:length(paras)
    [~,A1,A2] = make_nitrate_curve(paras(ii,:));
    a1s(ii,:) = A1;
    a2s(ii,:) = A2;
end
Dynamics_Fitted = struct('Time',[],'Nitrate_none',[],'Nitrate_chl',[]);
Dynamics_Fitted.Time = T;
Dynamics_Fitted.Nitrate_none = a1s;
Dynamics_Fitted.Nitrate_chl = a2s;

save('results/Fitted_curve.mat','Dynamics_Fitted');



function [T,A1,A2] = make_nitrate_curve(paras)

a0c = paras(1);
a0n = paras(2);
x0 = paras(3);
ga = 4.8;
ts = paras(4);

T = 0:0.001:4;

A1 = a0n-x0/ga.*(exp(ga*min(T,ts))-1)-x0*exp(ga*ts)*max(T-ts,0);
A1 = max(A1,0);
A2 = max(a0c-x0*T,0);


end