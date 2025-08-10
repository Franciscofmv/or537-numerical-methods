% 3rd order differentiation of a continuous function
function central_difference_example
f=@test_function;
x=-pi:0.5:pi; % differentiation at x=1, x is a vector
h=1e-3; %spacing or perturbation 10^(-2) to the -4 is recommended h.
n=3; % order of differentiation, we want to take third derivative
% uses central difference function (stored as separate file)
diffn=central_difference(f,x,h,n); % vector central differentiation
% diffn=richardson_diff(f,x,h,n);
plot(x,diffn,'-r'); % plotting numerical
hold on;
% analytical differentiation:
syms x;
diffa=double(subs(diff(x^2*exp(x),3),x,-pi:0.5:pi)); 
x=-pi:0.5:pi;
plot(x,diffa,'go'); % why we have go
legend('numerical','anaytical');
xlabel('x');
ylabel('d^3f/dx^3');
return

% test function:
function [f] = test_function(x)
f = x.^2.*exp(x);
return