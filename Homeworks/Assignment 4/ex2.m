function ex2
C=[5; 3; 0]; % initial condition
t= [0:0.012:1.2]; % 100 space vector
f=@odesyst; 
% 
% % Runge Kutta:
% [result]=runge_kutta(f,t,C);
% c1=result(1,:);
% c2=result(2,:);
% c3=result(3,:);
% 
% figure(1)
% plot(t,c1,'ro',t,c2,'-b+',t,c3,'-k*')
% title('Runge Kutta');
% legend('C1','C2','C3');
% xlabel('Time (minutes)')
% ylabel('Molar concentration  (moles/m^3)')
% 
% %Ode 45:
% figure(2)
% ode45(f,t,C)
% title('ode 45')
% legend('C1','C2','C3');
% xlabel('Time (minutes)')
% ylabel('Molar concentration  (moles/m^3)')
% 
%  
% figure(3)
% plot(t,c1,'-ro',t,c2,'-b+',t,c3,'-g^');
% xlabel('Time (minutes)')
% ylabel('Molar concentration  (moles/m^3)')
% hold on
% ode45(f,t,C)
% title('ode 45 and runge kutta');
% 
% % for new alpha's
% 
f2=@odesyst2;
% [result2]=runge_kutta(f2,t,C);
% c1=result2(1,:);
% c2=result2(2,:);
% c3=result2(3,:);
% 
% figure(4)
% plot(t,c1,'ro',t,c2,'-b+',t,c3,'-k*')
% title('Runge Kutta');
% legend('C1','C2','C3');
% xlabel('Time (minutes)')
% ylabel('Molar concentration  (moles/m^3)')
% xlim([0 1.2])
% %Ode 45:
% figure(5)
% ode45(f2,t,C)
% title('ode 45')
% legend('C1','C2','C3');
% xlabel('Time (minutes)')
% ylabel('Molar concentration  (moles/m^3)')
% xlim([0 1.2])

 
% figure(6)
% plot(t,c1,'-ro',t,c2,'-b+',t,c3,'-g^');
% xlabel('Time (minutes)')
% ylabel('Molar concentration  (moles/m^3)')
% xlim([0 1.2])
% hold on
ode23s(f2,t,C)
title('ode 45 and runge kutta');

end




%chemical_reactor function
function [dcdt] = odesyst(t,c)
c01=3.2;c02=4.8;Q=10; V=2;
dcdt=zeros(3,1);
dcdt(1)=(Q/V)*(c01-c(1))-2.6*c(1)*c(2);
dcdt(2)=(Q/V)*(c02-c(2))-2.6*c(1)*c(2);
dcdt(3)=(Q/V)*c(3)+2.6*c(1)*c(2);
end
function [dcdt] = odesyst2(t,c)
c01=3.2;c02=4.8;Q=10; V=2;
dcdt=zeros(3,1);
dcdt(1)=(Q/V)*(c01-c(1))-100000*c(1)*c(2);
dcdt(2)=(Q/V)*(c02-c(2))-1*c(1)*c(2);
dcdt(3)=(Q/V)*c(3)+10*c(1)*c(2);
end


function [y] = runge_kutta(df,x,y0)
m = length(y0);
n = length(x);
% y is a matrix of size 'mxn'
y=zeros(m,n);
y(:,1)=y0;
for i=1:n-1
    h = x(i+1)-x(i);
    xi=x(i);
    yi=y(:,i);
    k1 = feval(df,xi,yi);
    k2 = feval(df,xi+0.5*h,yi+0.5*k1*h);
    k3 = feval(df,xi+0.5*h,yi+0.5*k2*h);
    k4 = feval(df,xi+h,yi+k3*h);
    y(:,i+1) = y(:,i) + h*(k1+2*k2+2*k3+k4)/6;
end
end
