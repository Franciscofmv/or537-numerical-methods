function test_script
x=rand(100,1); % x is a column vector
[mean, stdev]=stat(x);
fprintf('mean = %g\n',mean);
fprintf('stdev = %g\n',stdev);
end

function [mean,stdev] = stat(x)
%STAT Interesting statistics.
n = length(x);
mean = avg(x,n);
stdev = sqrt(sum((x-avg(x,n)).^2)/n);
end

%-------------------------
function mean = avg(x,n)
%MEAN subfunction
mean = sum(x)/n;
end
