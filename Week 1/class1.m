% MATLAB can be used interactively
% file are format program.m 
% help window is available
a=1:2 % row vector
a=rand(5,1) % column vector of 5 random numbers between 0 and 1
a(4,1)
a
a'
a.^2 % since a is a vector
b=rand(5,1);
a.*b % element by element multiplication
a~=b
a = [1,2,3]% row vector
a = [1 2 3]% row vecotr
a=[1;2;3]% semicolon goes to the same row
% built in functions in MATLAB
zeros(10,1)
ones(5,1)
zeros(1,10)
ones(1,10)
a = 0:2:10 % row vector from 0 to 10 in increment of 2
zeros(5) % 5by5 matrix
a = rand(5,1)
A = zeros(5)
size(A) % 5 rows, and columns
a = 0:2:10
length(A)

eye(2)% identity matrix 
A=diag(1:3)
A = rand(5)
A(:,2)
max(A) % largest element of all columns
max(max(A))
%% Example 1%%
Q = [       5 4 6;
            3 2 4;
            6 5 3;
            3 5 4;
            2 4 3] % quantity in tons
p = [300;550;400;250;500] % price per ton
materials = 1:5 % 1 through 5 vector

% Total spent in May, June and July
%p = p' % now p is equal to p'
%total_spent_by_month = p*Q
total_spent_per_material = p.*sum(Q,2)
Total_spent = sum(total_spent_per_material)





