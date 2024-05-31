% analysis on class level, find Clostridia boundaries

clear
load('results/abundance_coarse_grain.mat');
load('results/function_and_phylum_data.mat','phase')

lv = 4;
nn = 1:118;
a00 = Abd_data(lv).Native(:,nn);
a2 = Abd_data(lv).None(:,nn);
a1 = Abd_data(lv).CHL(:,nn);
a0 = zeros(130,length(nn));
for ii=1:10
    kk = (1:13)+(ii-1)*13;
    a0(kk,:) = ones(13,1)*a00(ii,:);
end

lg0 = log(sum(a00,1)+1e-3)/log(10);
lg1 = log(a1+1e-3)/log(10)-log(a0+1e-3)/log(10);
lg2 = log(a2+1e-3)/log(10)-log(a1+1e-3)/log(10);
ph = reshape(ph1,[13,10]);
c24 = reshape(lg2(:,24),[13,10]);
phase = reshape(phase,[13,10]);
ys = 0.5;
bb = (c24>ys);

boundary_class_24 = NaN(10,1);
for ii=1:10
    for jj=1:12
        if (bb(jj,ii)==0)&&(bb(jj+1,ii)==1)
            a = (c24(jj+1,ii)-ys)/(c24(jj+1,ii)-c24(jj,ii));
            bc = ph(jj,ii)*a+ph(jj+1,ii)*(1-a);
            boundary_class_24(ii) = bc;
        end
    end
end

save('results/class_boundary.mat','c24','boundary_class_24');