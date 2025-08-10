syms x % create symbolic variables and functions
int(x) % integrate x, we will get a symbolic answer
diff(x) % differentiate x and get a symbolic answer
syms x c% define two variables 
int(c*x,x)% integrate c*x with respect to x
%%creating a symbolic ode equation%%
syms c(t) mu X K% creating symbolic function c(t)
ode = diff(c,t)== -4*0.5*(c/(12+c)) % define the ode
cond1 = c(0)==2; % 
sol1(t) = dsolve(ode,cond1)
