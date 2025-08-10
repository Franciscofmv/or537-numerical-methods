function int_bungee
syms g c m v t T; % defining symbolic variables
f(v) = sqrt(g*m/c)*tanh(sqrt(g*c/m)*t)% symbolic function
int_result = int(f,t,0,T) % integral with respect to t from 0 to T
g = 9.8; c = 0.1; m = 78; T = 10;
double(subs(int_result)) % this is my final result

end
