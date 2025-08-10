% chapra and canale p. 218
function beam_example
x0=[6;];
f=@displacement;
x_fsolve=fsolve(f,x0)
x_nr=newton_raphson(f,x0)
fplot(f,[0 20]);
grid on;
return


function [f] = shear(x)
f=20*(singularity(x,0)^1+singularity(x,5)^1)-15*singularity(x,8)^1-57;
return

function [f] = displacement(x)
f=-(5/6)*(singularity(x,0)^4-singularity(x,5)^4)+...
    (15/6)*singularity(x,8)^3+...
    75*singularity(x,7)^2+...
    (57/6)*x^3-238.25*x;
return

function [f] = singularity(x,a)
if (x > a)
    f=x-a;
else
    f=0;
end
return
