% fit nitrate data for model parameters

clear
load('results/model_parameters_fitting_nitrate.mat');
load('results/nitrate_data_for_fit_cleaned.mat');

options = optimoptions('fmincon','Display','off');
para0 = paras;
paras = para0;
sample_size = size(paras,1);
er = zeros(sample_size,1);
for ii=1:sample_size
    e1 = @(p)er_function(p,fdata{ii,1},fdata{ii,2});
    lb = [1.5;1.5;0;0];
    ub = [3.5;3.5;1;4];
    para = fmincon(e1,para0(ii,:),[],[],[],[],lb,ub,[],options);
    er(ii) = er_function(para,fdata{ii,1},fdata{ii,2});
    paras(ii,:) = para;
end

ts = paras(:,4);
x0 = paras(:,3);
a0 = paras(:,2);
ga = 4.8;
ts = min(log(a0*ga./x0+1),ts);
c0 = min(x0/ga.*(exp(ga*ts)-1),a0);
paras(:,4) = ts;
model_parameters(:,'error') = array2table(er);
model_parameters(:,'A0c') = array2table(paras(:,1));
model_parameters(:,'A0n') = array2table(paras(:,2));
model_parameters(:,'x0') = array2table(x0);
model_parameters(:,'C0') = array2table(c0);
model_parameters(:,'t_star') = array2table(ts);

save('results/model_parameters_fitting_nitrate.mat','model_parameters','paras');

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

er1 = mean((a1-A1).^2);
er2 = mean((a2-A2).^2);
er = (er1+er2)/(length(t1)+length(t2));

end