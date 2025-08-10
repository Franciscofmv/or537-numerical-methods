% palm 2_64 example
location = {'Smith. St.','Hope Ave', 'Clark St.', 'North Rd'};
max_load = {80, 90, 85, 100};
year_built = {1928, 1950, 1933, 1960};
due_maintenance = {1997, 1999, 1998, 1998};
bridge = struct('location',location,'max_load',max_load,...
        'year_built',year_built,'due_maintenance',due_maintenance);
getfield(bridge,{2},'location')
names = fieldnames(bridge)
% palm 2_65
bridge(3).due_maintenance=2000;
getfield(bridge,{3},'due_maintenance')
% palm 2_66
bridge(5).location = 'Shore Rd.';
bridge(5).max_load = 85;
bridge(5).year_built = 1997;
bridge(5).due_maintenance = 2002;
bridge
%end


