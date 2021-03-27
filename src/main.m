clc;
clear;
close all;

%% Tests sur le premier type de maillage constitué uniquement d'éléments triangles + conditions de Dirichlet uniquement
U1 = elliptic(5, @f, @g, false, true);

%% Inclusion du traitement de tels éléments de type Q1 et des conditions de Neumann
U2 = elliptic(0, @f, @g, true, true);

%% l'évolution en loi log-log de || u_ex - u_h ||_2h en focntion de h 
n_max = 30;
y = zeros(n_max - 9, 1);
x = zeros(n_max - 9, 1);

for i = 10:n_max
    uh  = elliptic(i, @f2, @g, false, false);
    h    = 1 / sqrt(length(uh));
    [coord, ~] = maillage_carre(i);
    uex = u_ex(coord(:, 1), coord(:, 2));
    x(i - 9) = h;
    y(i - 9) = h * norm(uex - uh, 2);
end

figure(3)
loglog(x, y);
grid on
ylabel ('|| u_ex - u_h ||_2h');
xlabel('h');
title('évolution en loi log-log de || u_{ex} - u_h ||_{2h} en focntion de h ');

%% l'évolution du nombre d'éléments non nuls de R en fonctionde la taille de la matrice
