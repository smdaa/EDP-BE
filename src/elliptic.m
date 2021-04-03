function [U, A] = elliptic(n, f, g, Q1, aff)  
%   n   = nombre de points par cote du carré
%   f   = fonction f
%   g   = fonction g
%   Q1  = inclusion du traitement des éléments de type Q1 et des conditions de Neumann 
%   aff = affichage de la solution

    if Q1
        elements4 = load('elements4.dat');
        elements3 = load('elements3.dat');
        neumann = load('neumann.dat');
        coordinates = load('coordinates.dat');
        dirichlet = load('dirichlet.dat');
    else
        [coordinates, elements3, elements4, dirichlet, ~] = maillage_carre(n);
    end

    Ns = size(coordinates, 1);  % Nombre de points du maillage
    
    % Assemblage de la matrice A et du second membre b dans le cas d’un maillage constitue uniquement d'elements triangles
    A = zeros(Ns);
    B = zeros(Ns, 1);
    for ind = elements3'
        Tloc = [coordinates(ind(1), :); coordinates(ind(2), :); coordinates(ind(3), :)];
        [M, alpha] = raideur_triangle(Tloc);
        A(ind, ind) = A(ind, ind) + M;
        B(ind) = B(ind) + (alpha / 6) .* f((sum(coordinates(ind, :)) ./ 3));
    end

    % Contribution des conditions de Dirichlet
    for i=1:Ns
           B(i) = B(i) - A(i, dirichlet) * u_d(coordinates(dirichlet, :));
    end

    % Assemblage de la matrice A et du second membre b  dans le cas d’un maillage constitue d'elements Quadrangles;
    if Q1
        for ind = elements4'
            Qloc = [coordinates(ind(1),:); coordinates(ind(2),:); coordinates(ind(3),:); coordinates(ind(4),:)];
            [M, alpha] = raideur_quadrangle(Qloc);
            A(ind, ind) = A(ind, ind) + M;
            B(ind) = B(ind) + (alpha / 4) .* f((sum(coordinates(ind, :)) ./ 4));
        end

        % Contribution des conditions de Neumann
        for ind = neumann'  
           a1 = coordinates(ind(1), :);
           a2 = coordinates(ind(2), :);
           B(ind) = B(ind) + (norm(a2 - a1)) * g((a1 + a2) / 2); 
        end
    end

    U = zeros(Ns, 1);
    % condittion de Dirichlet
    U(dirichlet) = u_d(coordinates(dirichlet));
    temp = setdiff(1:Ns, dirichlet); 
    U(temp) = A(temp,temp) \ B(temp);
    if aff
        %%Affichage
        figure;
        show(elements3, elements4, coordinates, U);
    end
end