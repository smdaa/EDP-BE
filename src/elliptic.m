clc
n=3
[coordinates, elements3, elements4, dirichlet, neumann] = maillage_carre(n)

%Tloc=[coordinates(1,:);coordinates(2,:); coordinates(4,:)]
%Mloc=raideur_triangle(Tloc)

Ns=n ^ 2
A=zeros(Ns)

for ind=elements3'
    ind
    Tloc=[coordinates(ind(1),:);coordinates(ind(2),:); coordinates(ind(3),:)]
    M=raideur_triangle(Tloc)
    A(ind,ind)=A(ind,ind)+M
end

B = zeros(Ns,1);

for ind=elements3'
    T=[coordinates(ind(1),:);coordinates(ind(2),:); coordinates(ind(3),:)];
    alpha=det([ (T(2,1)-T(1,1)) , (T(3,1)-T(1,1)) ; (T(2,2)-T(1,2)) , (T(3,2)-T(1,2)) ]);
    B(ind) = B(ind) + (alpha / 6) .* f((sum(coordinates(ind, :)) ./ 3));
        
end
B
    

    
        

