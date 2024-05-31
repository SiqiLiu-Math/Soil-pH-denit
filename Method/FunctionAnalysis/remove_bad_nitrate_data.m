% Remove bad nitrate data manually

clear
load('results/nitrate_data_for_fit.mat');

% [sample, CHL-/+, timepoint]
b5n = [(157:195)',ones(39,1),ones(39,1)*2];
b5c = [(157:195)',ones(39,1)*2,ones(39,1)*2];
b7n = [(235:273)',ones(39,1),ones(39,1)*3];
b7c = [(235:273)',ones(39,1)*2,ones(39,1)*3];
b12n = [(430:468)',ones(39,1),ones(39,1)*3];
b12c = [(430:468)',ones(39,1)*2,ones(39,1)*3];
bad_data = [b5c;b5n;b7c;b7n;b12c;b12n]; % bad data in soil 5/7/12
otherbad = [338,2,4;204,1,8;382,1,4;714,2,2];
bad_data = [bad_data;otherbad];

for ii=1:size(bad_data,1)
    fd = fdata{bad_data(ii,1),bad_data(ii,2)};
    fd(:,bad_data(ii,3)) = [];
    fdata{bad_data(ii,1),bad_data(ii,2)} = fd;
end

save('results/nitrate_data_for_fit_cleaned.mat','fdata')