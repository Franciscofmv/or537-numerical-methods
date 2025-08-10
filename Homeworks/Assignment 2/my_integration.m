function my_integration
a = 0; % lower limit
b = 1; % upper limit
n=10
integrand(1)

end


% 
% 
% function [I] = composite_trapezoidal(f,a,b,n)
% I = 0;
% h = (b-a)/n;
%     for i = 0:n
%         x = a + i*h;
%         if (i==0 || i==n)
%             I = I + f(x)
%         else
%             I = I + 2*f(x)
%         end 
%     end
%     I = I*h/2
% end

function [y] = integrand(x)
y=x.^(0.1)*(1.2-x)*(1-exp(x))
end
