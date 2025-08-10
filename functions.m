function test_script
  x = rand(100,1); % random column vector
  [mean,stdev]=stat(x);
  fprintf('mean = %g\n',mean)
  fprintf('stdev= %g\n',stdev)
 end


function [mean,stdev] = stat(x)
  % this functions gives an 
  n = length(x);
  mean = avg(x,n); %calling another function inside this function
  stdev = sqrt(sum((x-avg(x,n)).^2)/n);
  end
  
  function mean = avg(x,n)
  % calculate the mean of the elements of a vector
  mean=sum(x)/n;
  end