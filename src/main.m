clc;
clear;
close all;

%% Tests sur le premier type de maillage constitué uniquement d'éléments triangles + conditions de Dirichlet uniquement
elliptic(5, @f, @g, false, true);

%% Inclusion du traitement de tels éléments de type Q1 et des conditions de Neumann
elliptic(0, @f, @g, true, true);

%% l'évolution en loi log-log de || u_ex - u_h ||_2h en focntion de h 
n_max = 30;
y  = zeros(n_max - 9, 1);
x  = zeros(n_max - 9, 1);

for i = 10:n_max
    [uh, ~]  = elliptic(i, @f2, @g, false, false);
    h    = 1 / sqrt(length(uh));
    [coord, ~] = maillage_carre(i);
    uex = u_ex(coord(:, 1), coord(:, 2));
    x(i - 9) = h;
    y(i - 9) = h * norm(uex - uh, 2);
end

y2 = x .^ 2;
y3 = x .^ 3;
y1 = x;

figure(3)
loglog(x, y);
hold on
loglog(x, y1);
hold on
loglog(x, y2);
hold on
loglog(x, y3);
grid on
ylabel ('|| u_ex - u_h ||_{2h}');
xlabel('h');
legend('|| u_ex - u_h ||_{2h}', 'h', 'h^2','h^3')
title('évolution en loi log-log de || u_{ex} - u_h ||_{2h} en focntion de h ');

%% l'évolution du nombre d'éléments non nuls de R en fonctionde la taille de la matrice
n_max = 30;

x = zeros(n_max - 9, 1);
y = zeros(n_max - 9, 1);

for i = 10:n_max
    [U, A]  = elliptic(i, @f, @g, false, false);
    [V, D] = eig(A);
    D(D < 0) = 0.0;
    A = V * D * (V'); 
    R = chol(abs(A));
    x(i - 9) = length(A);
    y(i - 9) = length(find(R));    
end

figure(4)
plot(x, y)
grid on
ylabel ('nombre d éléments non nuls de R');
xlabel('taille de la matrice A');
title('évolution du nombre éléments non nuls de R en fonctionde la taille de la matrice');