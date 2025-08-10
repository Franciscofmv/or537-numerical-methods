% palm 2.61 (p.118)
% cell array example
A = cell(2,2); % initialize cell array A
A{1,1} = 'Motor 28C'; 
A{1,2} = 'Test ID 6';
A{2,1} = [3 9; 7 2];
A{2,2} = [6 5 1];
celldisp(A); % dsiplay cell contents
cellplot(A); % plot cell array
%end
