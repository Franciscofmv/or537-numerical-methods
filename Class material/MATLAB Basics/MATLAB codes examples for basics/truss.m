% truss example C&C p.308
% define truss matrix
A=zeros(6); b=zeros(6,1);
A=[0.866 0 -0.5 0 0 0;
    0.5 0 0.866 0 0 0;
    -0.866 -1 0 -1 0 0;
    -0.5 0 0 0 -1 0;
    0 1 0.5 0 0 0;
    0 0 -0.866 0 0 -1;];
% rhs for original case
b=[0; -1000; 0; 0; 0; 0;];
x=round(A\b);
disp('ORIGINAL CASE');
disp('-------------');
fprintf('F1 = %g\n',x(1));
fprintf('F2 = %g\n',x(2));
fprintf('F3 = %g\n',x(3));
fprintf('H2 = %g\n',x(4));
fprintf('V2 = %g\n',x(5));
fprintf('V3 = %g\n\n',x(6));
% case 2
C=inv(A);
b=[1000; 0; 1000; 0; 0; 0;];
x=round(C*b);
disp('CASE 2');
disp('------');
fprintf('F1 = %g\n',x(1));
fprintf('F2 = %g\n',x(2));
fprintf('F3 = %g\n',x(3));
fprintf('H2 = %g\n',x(4));
fprintf('V2 = %g\n',x(5));
fprintf('V3 = %g\n\n',x(6));
% case 3
b=[-1000; 0; 0; -1000; 0; 0;];
x=round(C*b);
disp('CASE 3');
disp('------');
fprintf('F1 = %g\n',x(1));
fprintf('F2 = %g\n',x(2));
fprintf('F3 = %g\n',x(3));
fprintf('H2 = %g\n',x(4));
fprintf('V2 = %g\n',x(5));
fprintf('V3 = %g\n',x(6));
