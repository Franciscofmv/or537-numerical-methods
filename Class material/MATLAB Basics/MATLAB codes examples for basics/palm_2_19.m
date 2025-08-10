% cost estimation example - palm 2.19 (p.100)
% Q: 5 x 3 matrix of quantity purchased
Q = [5 4 6;
    3 2 4;
    6 5 3;
    3 5 4;
    2 4 3];
% p: 5 x 1 vector of price of each material
p = [300; 550; 400; 250; 500];
% transpose p into a row vector of 1 x 5
% and multiply by matrix Q
cost_per_month = p'*Q
% sum Q along each row and dot multiply by p
cost_per_material = p.*sum(Q,2)
% total spent
total = sum(cost_per_material)
