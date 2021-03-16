clc;
clear;
close all;
%% Initialisation
Q1 = true;      % Traitement Q1
n = 8;          % Nombre de points par cote du carré
if Q1
    elements4 = load('elements4.dat');
    elements3 = load('elements3.dat');
    neumann = load('neumann.dat');
    coordinates = load('coordinates.dat');
    dirichlet = load('dirichlet.dat');
else
    [coordinates, elements3, elements4, dirichlet, ~] = maillage_carre(n);
end

Ns = size(coordinates, 1);
%% Assemblage de la matrice A dans le cas d’un maillage constitue uniquement d'elements triangles;
A = zeros(Ns);
for ind = elements3'
    Tloc = [coordinates(ind(1), :); coordinates(ind(2), :); coordinates(ind(3), :)];
    [M, alpha] = raideur_triangle(Tloc);
    A(ind, ind) = A(ind, ind) + M;
end

%% Assemblage du second membre dans le cas de conditions de Dirichlet uniquement
B = zeros(Ns, 1);
for ind = elements3'
    B(ind) = B(ind) + (alpha / 6) .* f((sum(coordinates(ind, :)) ./ 3));     
end

% Contribution des conditions de Dirichlet
for i=1:Ns
       B(i) = B(i) - A(i, dirichlet) * u_d(coordinates(dirichlet, :));
end

%% Assemblage de la matrice A dans le cas d’un maillage constitue d'elements Quadrangles;
if Q1
    for ind = elements4'
        Qloc = [coordinates(ind(1),:); coordinates(ind(2),:); coordinates(ind(3),:); coordinates(ind(4),:)];
        M = raideur_quadrangle(Qloc);
        A(ind, ind) = A(ind, ind) + M;
    end

%% Assemblage du second membre dans le cas de conditions de Dirichlet et de Neumann
    for ind = elements4'
        B(ind) = B(ind) + (alpha / 4) .* f((sum(coordinates(ind, :)) ./ 4));     
    end

    % Contribution des conditions de Neumann
    for ind = neumann'  
       a1 = coordinates(ind(1), :);
       a2 = coordinates(ind(2), :);
       B(ind) = B(ind) + (norm(a2 - a1) / 2) * g((a1 + a2) / 2); 
    end
end

U = zeros(Ns, 1);
% condittion de Dirichlet
U(dirichlet) = u_d(coordinates(dirichlet));
temp = setdiff(1:Ns, dirichlet); 
U(temp) = A(temp,temp) \ B(temp);
figure;
show(elements3, elements4, coordinates, U);
