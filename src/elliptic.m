clc
n = 5;
[coordinates, elements3, elements4, dirichlet, neumann] = maillage_carre(n);

Ns = n ^ 2;

%Assemblage de la matrice A
A = zeros(Ns);
for ind = elements3'
    Tloc = [coordinates(ind(1),:);coordinates(ind(2),:); coordinates(ind(3),:)];
    [M,alpha] = raideur_triangle(Tloc);
    A(ind,ind) = A(ind,ind) + M;
end

%Assemblage du vecteur B
B = zeros(Ns,1);
for ind = elements3'
    T = [coordinates(ind(1),:);coordinates(ind(2),:); coordinates(ind(3),:)];
    B(ind) = B(ind) + (alpha / 6) .* f((sum(coordinates(ind, :)) ./ 3));     
end

% Contribution des conditions de Dirichlet
for i=1:Ns
       B(i) = B(i) - A(i,dirichlet) * u_d(coordinates(dirichlet,:));
end

U = zeros(Ns, 1);
U(dirichlet) = u_d(coordinates(dirichlet));
temp = setdiff(1:Ns, dirichlet);   
U(temp) = A(temp,temp) \ B(temp);
show(elements3, elements4, coordinates, U);

