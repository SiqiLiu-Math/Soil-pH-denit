% Find the turning point of phylum2 and phylum9 death and phlyum4 growth

clear
load('results/phylum_pca_nmf_data.mat');
ph = reshape(ph1,[13,10]);
ph0b = ph0(1:13:end);
p2 = reshape(lg1(:,2),[13,10]);
ys2 = -0.2;
l2 = zeros(13,10);
p9 = reshape(lg1(:,9),[13,10]);
ys9 = -0.4;
l9 = zeros(13,10);
p4 = reshape(lg2(:,4),[13,10]);
ys4 = 0.5;
l4 = zeros(13,10);

for ii=1:10
    l2(p2(:,ii)>ys2,ii) = 1;
    l9(p9(:,ii)>ys9,ii) = 1;
    l4(p4(:,ii)>ys4,ii) = 1;
end
bd = zeros(10,3);
for ii=1:10
    for jj=1:12
        s2 = (l2(jj,ii)==0)+(l2(jj+1,ii)==1);
        s9 = (l9(jj,ii)==0)+(l9(jj+1,ii)==1);
        s4 = (l4(jj,ii)==0)+(l4(jj+1,ii)==1);
        if s2==2
            a = (p2(jj+1,ii)-ys2)/(p2(jj+1,ii)-p2(jj,ii));
            bd(ii,1) = ph(jj,ii)*a+ph(jj+1,ii)*(1-a);
        end
        if s9==2
            a = (p9(jj+1,ii)-ys9)/(p9(jj+1,ii)-p9(jj,ii));
            bd(ii,2) = ph(jj,ii)*a+ph(jj+1,ii)*(1-a);
        end
        if s4==2
            a = (p4(jj+1,ii)-ys4)/(p4(jj+1,ii)-p4(jj,ii));
            bd(ii,3) = ph(jj,ii)*a+ph(jj+1,ii)*(1-a);
        end
    end
end

load('../Function_fitting_double_cr_model/results/function_parameter_and_phases.mat');
p1 = [1,2,3,2,3,4,3,3,3,2,3,5,7,8,8,9,9,10,5];
for ii=1:19
    jj = (ii-1)*13+(1:p1(ii));
    phase(jj)=1;
end
ll = length(x0);
fd = zeros(20,8);
unt = -unt;
for ii=1:ll-1
    s12 = (phase(ii)~=2)+(phase(ii+1)==2);
    s23 = (phase(ii)==2)+(phase(ii+1)==3);
    if s12==2
        nn = sid(ii+1);
        fd(nn,1) = ph1(ii+1)*0.5+min(ph1(ii+1),ph1(ii))*0.5;
        fd(nn,3) = (ph1(ii+1)-min(ph1(ii+1),ph1(ii)))*0.5;
        fd(nn,5) = unt(ii+1)*0.5+min(unt(ii+1),unt(ii))*0.5;
        fd(nn,7) = (unt(ii+1)-min(unt(ii+1),unt(ii)))*0.5;
    end
    if s23==2
        fd(nn,2) = ph1(ii+1)*0.5+ph1(ii)*0.5;
        fd(nn,4) = (ph1(ii+1)-ph1(ii))*0.5;
        fd(nn,6) = unt(ii+1)*0.5+unt(ii)*0.5;
        fd(nn,8) = (unt(ii+1)-unt(ii))*0.5;
    end
end
[~,uu] = unique(sid);
ph0f = ph0(uu);

x1 = [ones(length(ph0f),1),ph0f];
b1 = x1(1:end-1,:)\fd(1:end-1,1);
b2 = x1(1:end-1,:)\fd(1:end-1,2);
x2 = [ones(length(ph0b),1),ph0b];
b3 = [x2;x2]\[bd(:,1);bd(:,2)];
b4 = x2\bd(:,3);
b5 = [x1;x2;x2]\[fd(:,1);bd(:,1);bd(:,2)];
b6 = [x1;x2]\[fd(:,2);bd(:,3)];
line_fit = struct('Function12',b1,'Function23',b2,'Abundance12',b3,'Abundance23',b4,'All12',b5,'All23',b6);

save('results/phase_boundaries2.mat','fd','bd','ph0b','ph0f','line_fit');