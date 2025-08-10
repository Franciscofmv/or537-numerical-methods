% palm 2_62 cell array example
A = cell(2,4);
A{1,1} = 'distance';
A{1,2} = 'Length';
A{1,3} = 'radius';
A{1,4} = 'Capacitance';
d=[0.003,0.004,0.005,0.01];
L=[1,2,3];
r=[0.001,0.002, 0.003];
A{2,1} = d;
A{2,2} = L;
A{2,3} = r;
eps = 8.854e-12;
C = zeros(4,3,3);
for i=1:4
    for j=1:3
        for k=1:3
            C(i,j,k)=pi*eps*L(j)/(log((d(i)-r(k))/r(k)));
        end
    end
end
A{2,4} = C;
celldisp(A);
A{2,4}(3,2,1)
cellplot(A);
% end