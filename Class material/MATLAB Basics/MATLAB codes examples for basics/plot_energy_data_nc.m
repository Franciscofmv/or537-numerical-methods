% plot energy data for NC
function plot_energy_data_nc
A=readcell('energy_generation_nc.xlsx','Range','B6:C10');
nrows = size(A,1);
labels = cell(nrows,1);
values = zeros(nrows,1);
for i=1:nrows
    labels{i} = A{i,2};
    values(i) = A{i,1};
end
nc_energy = struct('lables',labels,'values',values);
save('ncenergy.mat','nc_energy')
% plot results
figure(1);
pie(values);
legend(labels,'Location','eastoutside','Orientation','vertical');
title(['Energy generation by source in 2014 for NC']);
end