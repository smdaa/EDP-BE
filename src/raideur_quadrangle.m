function [M, alpha] = raideur_quadrangle(Q)
    alpha = det([(Q(2,1) - Q(1,1)) , (Q(3,1) - Q(1,1)) ; (Q(2,2) - Q(1,2)) , (Q(3,2) - Q(1,2))]);
    Jac_phi = [Q(2, 1) - Q(1, 1), Q(4, 1) - Q(1, 1); Q(2, 2) - Q(1, 2), Q(4, 2) - Q(1, 2)];
    temp = pinv(Jac_phi' * Jac_phi);
    
    a = temp(1, 1);
    b = temp(1, 2);
    c = temp(2, 1);
        
    M = (2 * (a + c) + 3 * b) * eye(4);
    M(1:3, 4) = [-2 * c + a; -(a + 3 * b + c); c - 2 * a];
    M(1:2, 3) = [-(a + 3 * b + c); -2 * c + a];
    M(1, 2) = -2 * c + a;
    M = floor((M + M.') / 2);
end
