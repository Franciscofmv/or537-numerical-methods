% tensile force example palm 2_24 (p. 102)
% define scalar constants
Lb = 5; Lc = 5; W = 400;
% define D vector
D = 0.1:0.1:4.9;
% calculate T vector
T = Lb*Lc*W./(D.*sqrt(Lc^2 - D.^2));
% plot T vs D
plot(D,T);
% find minium of T and the corresponding index
[min_T, min_index] = min(T);
% print the D value corresponding to minium T
D_min = D(min_index)
% find values of D within 10% of min
[temp, design_index] = find(T < 1.1*min_T);
% extract D values that correspond to T_design
D_design = D(design_index);
T_design = T(design_index);
plot(D_design,T_design);
% end
