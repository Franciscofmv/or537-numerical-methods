% chapra and canale 21.12 (p.541)
function jet_fighter
t = [0 0.52 1.04 1.75 2.37 3.25 3.83];
x = [153 185 210 249 261 271 273];
v = tabular_derivative(t,x);
figure(1);
plot(t,v);
xlabel('time (s)');
ylabel('velocity (m/s)');
a = tabular_derivative(t,v);
figure(2);
plot(t,a);
xlabel('time (s)');
ylabel('acceleration (m/s^2)');
return
