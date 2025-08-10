syms f(c) k V F Q;
f(c)=c^2*k*V+Q*c-F
V=12;F=175;Q=1;k=0.15;
fun=subs(f)
sol = solve(fun,c)
final=double(sol(2));
%%
 syms g(c) k V F Q;
 
g(c)=V/(F-Q*c-k*V*c^2);
V=12;F=175;Q=1;k=0.15;c=final;
fun1 = subs(g)
double(fun1)

