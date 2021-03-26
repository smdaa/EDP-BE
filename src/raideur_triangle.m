function [M, alpha] = raideur_triangle(T)
    alpha = det([(T(2,1) - T(1,1)) , (T(3,1) - T(1,1)) ; (T(2,2) - T(1,2)) , (T(3,2) - T(1,2))]);
    M = zeros(3,3);
    for i = 0:2 
        for j = 0:2
            grad_nj = (1/alpha) * [(T(mod(j+1, 3)+1, 2) - T(mod(j+2, 3)+1, 2)); (T(mod(j+2, 3)+1, 1) - T(mod(j+1, 3) + 1, 1))];
            grad_ni = (1/alpha) * [(T(mod(i+1 ,3)+1, 2) - T(mod(i+2, 3)+1, 2)); (T(mod(i+2, 3)+1, 1) - T(mod(i+1, 3) + 1, 1))];
            M(i+1,j+1) = (alpha / 2) * grad_ni' * grad_nj;
        end
    end
end
