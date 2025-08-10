function [] = io_examples
n=5;
a=rand(1,n);
b=rand(1,n);
c=zeros(1,2*n);
% by default file is opened in matlab/work directory
file_id=fopen('random_vectors','w');
fwrite(file_id, a, 'double');
fwrite(file_id, b, 'double');
fclose(file_id);
fid=fopen('random_vectors','r');
[c, count]=fread(fid, 2*n, 'double');
fclose(fid);
sum_a=sum(a);
sum_b=sum(b);
sum_c=sum(c);
fprintf('sum of a %g\n',sum_a);
fprintf('sum of b %g\n',sum_b);
fprintf('sum of a + b %g\n',sum_a+sum_b);
fprintf('sum of c %g\n',sum_c);
return
