% basic nonlinear backward euler
%varargin = variable argument in, you can supply many arguments
function [y] = nonlinear_backward_euler_vararg(varargin)
if (nargin < 3)
    fprintf('need to supply atleast F, t, y0\n');
    return
else
    Ffunc = varargin{1}; % to first argument
    t = varargin{2}; % to second argument
    y0 = varargin{3};% 
end
n=length(t);
m=length(y0);
y=zeros(m,n);
y(:,1)=y0;
I=eye(m);
if (nargin > 3) % if n = 3 , user have not supplied jacobian
    Jfunc = varargin{4};
    for i=1:n-1
        dt=t(i+1)-t(i);
        J=feval(Jfunc,t(i+1),y(:,i));
        F=feval(Ffunc,t(i+1),y(:,i));
        dy=(I-dt*J)\(dt*F);
        y(:,i+1)=y(:,i)+dy;
    end
else
    Jfunc = @numerical_jacobian; % calculated using numerical jacobian func,separated code
    for i=1:n-1
        dt=t(i+1)-t(i);
        J=Jfunc(Ffunc,t(i+1),y(:,i));
        F=ffunc(t(i+1),y(:,i));
        dy=(I-dt*J)\(dt*F);
        y(:,i+1)=y(:,i)+dy;
    end
end
return
