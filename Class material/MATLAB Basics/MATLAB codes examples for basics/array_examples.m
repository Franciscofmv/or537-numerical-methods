function array_examples
% vector operations
vector_example;
matrix_example;
multidim_array_example;
cell_array_example;
structure_array_example;
return

function vector_example
%% row vectors %%%
x=[1 2 3 4 5] %create a 1x5 row vector x
x=[1,2,3,4,5] 
x=1:5                
x(1:5)=1:5          
% initialization
x=zeros(1,5)
x=ones(1,5)
x=rand(1,5)
% column vectors %
y=[1; 2; 3; 4; 5] % create a 5x1 columns vector y
y=[1
    2
    3
    4
    5]
y=x' % create column vector by transposing row vector
y=zeros(5,1)
y=ones(5,1)
y=rand(5,1)
z=x*y
z=x.*x
z=x*x'
z=dot(x,x)
%a row vector with complex numbers
x= [1+5*i 2+i 3+2*i]  % 1x3 vector
y=conj(x)
z=x.*y % multiply by its complex conjugate
return

% 2-dimensional matrices
function matrix_example
%create a 2x3 matrix x
x=[1 2 3; 4 5 6;]
x=[1 2 3
    4 5 6]
x=[1,2,3
    4,5,6]
x=[1, 2, 3; 4, 5, 6]
x=[1:3; 4:6]
y=x.^2 % result is a 3x2 matrix
y=x.*x  % same as above
y=x/2    % divide by a scalar
y=y'     % transpose x into 2x3 matrix
x=zeros(2,3)
x=rand(2,3)
x=ones(2,3)
y=x(1,:)  % store first row of x into y (3x1 vector)
z=x(:,2)  % store second column of x into z (1x2 vector)
return

% multidimensional arrays
function multidim_array_example
a=zeros(5,3,2) % 5x3x2 zero array
a=ones(5,3,2)
a=rand(5,3,2)
z=a(2,:,:)  % 3x2 matrix
z=a(1,:,1)  % 3x1 vector
b=a.*a      % 5x3x2 array
return

% cell arrays
function cell_array_example
A(1,1) = {'Neuse river basin'};
A(1,2) = {'Jan 10, 2005'};
A(2,1) = {[50, 60, 58]};
A(2,2) = {[45, 47, 46; 44, 46, 45; 42, 45, 43]};
celldisp(A);
cellplot(A);
return

function structure_array_example
student(1).name = 'John Smith';
student(1).SSN = '392-77-1786';
student(1).email = 'smithj@ncsu.edu';
student(1).tests = [67,75,84];
% display first elements
student
student(2).name = 'Mary Jones';
student(2).SSN = '431-56-9832';
student(2).email = 'jonesm@ncsu.edu';
student(2).tests = [84,78,93];
student(2)
% another way of initializing structures
name = {'John Smith', 'Mary Jones'};
email = {'smithj@ncsu.edu','jonesm@ncsu.edu'};
stud = struct('name',name,'email',email);
stud(1)
stud(2)
return
