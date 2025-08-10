function tab_ex3
T = [0 1 5.5 10 12 14 16 18 20 24]; % time in hours
C = [1 1.5 2.3 2.1 4 5 5.5 5 3 1.2]; % outflow concentration (mg/L)
flow = zeros(1,length(T));
    for i=1:length(T)
        flow(i) = flow_rate(i);
    end
 QC = flow.*C;
 g = @trap_tabular
 numerator = g(T,QC);
 denominator = integral(@flow_rate,0,24);
 Average_Concentration = numerator/denominator
end

function [Q] = flow_rate(t) % flow rate equation m^3/s
Q = 20 + 10* sin((t-10)*2*pi/24);
end