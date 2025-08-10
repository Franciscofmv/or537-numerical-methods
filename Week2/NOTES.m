name = {'John Smith','Mary Jones'}; % create cell arrays
SSN = {'392-77-1786','431-56-9832'};
email = {' smithj@ncsu.edu','jonesm@ncsu.edu'};
student = struct('name',name,'SSN',SSN,'email',email) % use command struct
 % access a structure field:
 student(1).name %struct_name(number of row-field).name of field
 student(2).SSN