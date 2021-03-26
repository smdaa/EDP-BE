clc;
clear;
close all;

n_max = 20;
y = zeros(n_max - 1, 1);
x = zeros(n_max - 1, 1);

for i = 2:n_max
    [coord, ~] = maillage_carre(i);
    uh = elliptic(false, i, false);
    h2 = 1/length(uh);
    x(i) = sqrt(h2);
    y(i) = h2 * norm(u_ex(coord(:, 1), coord(:, 2)) - uh, 2);
    
end
loglog(x, y);