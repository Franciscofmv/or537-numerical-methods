function [] = flow_control_examples
% some examples of for, if, while, and switch
for_example;
if_example;
while_example;
break_example;
switch_example;
return

function [] = for_example
disp('for_example');
disp('-----------');
n=10; m=5;
a=rand(n,m); b=rand(m,n); c=zeros(n,n);
% perform matrix product by hand
for i=1:n
    for j=1:n
        % dot product of ith row of a and jth column of b
        c(i,j)=c(i,j)+a(i,:)*b(:,j);
    end
end
% compute product
norm_calc=norm(c);
fprintf('norm of matrix product calculated by hand = %g\n',norm_calc);
% perform matrix product using "*" and compute norm for check
norm_matlab=norm(a*b);
fprintf('norm of matrix product calculated by simple multiplication = %g\n\n',norm_matlab);
return

function [] = if_example
disp('if_example');
disp('----------');
n=10000; p=0; q=0; r=0; s=0; t=0;
% generate random normal distribution between [-inf, inf]
% with mean zero and variance one into a row vector of 1xn
a=randn(1,n);
% do some arbitrary meaningless stuff
for i=1:n
    if (a(i) < -2 | a(i) > 2)
        p=p+1;
    elseif ((a(i) > -2 & a(i) < -1) | (a(i) > 1 & a(i) < 2))
        q=q+1;
    elseif (a(i) == -1 | a(i) == 0 | a(i) == 1)
        r=r+1;
    elseif (a(i) > -1 & a(i) < 1)
        s=s+1;
    else
        t=t+1;
    end
end
fprintf('percentage of numbers lower than -2 or greater than 2 = %g %%\n',100*p/n);
fprintf('percentage of numbers in ranges [-2,-1] and [1,2] = %g %%\n',100*q/n);
fprintf('percentage of numbers exactly equal to -1 or 0 or 1 = %g %%\n',100*r/n);
fprintf('percentage of numbers between -1 and 1 = %g %%\n',100*s/n);
fprintf('percentage of remaining numbers = %g %% \n\n',100*t/n);
return

function [] = while_example
disp('while_example');
disp('-------------');
% calculate pi using montecarlo method
in=0; out=0; k=0; eps=100;
% limit number of iterations to 10000 and tolerance to 10-4
while (k < 10000 & eps > 1e-4)
    x=2*(rand-0.5);
    y=2*(rand-0.5);
    r=sqrt(x^2+y^2);
    if (r < 1)
        in=in+1;
    else
        out=out+1;
    end
    k=k+1;
    pi_calc=4*in/(in+out);
    eps=abs(pi-pi_calc);
end
fprintf('number of iterations = %g\n',k);
fprintf('calculated pi = %g\n\n',pi_calc);
return


function [] = break_example
disp('break_example');
disp('-------------');
% calculate pi using montecarlo method
in=0; out=0;
% limit number of iterations to 10000 and tolerance to 10-4
for k=1:100000
    x=2*(rand-0.5);
    y=2*(rand-0.5);
    r=sqrt(x^2+y^2);
    if (r < 1)
        in=in+1;
    else
        out=out+1;
    end
    pi_calc=4*in/(in+out);
    eps=abs(pi-pi_calc);
    if (eps < 1e-4), break, end
end
fprintf('number of iterations = %g\n',k);
fprintf('calculated pi = %g\n\n',pi_calc);
return

function [] = switch_example
method = 'cubic';
switch method
  case {'linear','bilinear'}
    disp('Method is linear')
  case 'cubic'
    disp('Method is cubic')
  case 'nearest'
    disp('Method is nearest')
  otherwise
    disp('Unknown method.')
end
return

