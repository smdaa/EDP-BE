function [M, alpha] = raideur_quadrangle(Q)
    alpha = det([(Q(2,1) - Q(1,1)) , (Q(3,1) - Q(1,1)) ; (Q(2,2) - Q(1,2)) , (Q(3,2) - Q(1,2))]);
    Jac_phi = [Q(2, 1) - Q(1, 1), Q(4, 1) - Q(1, 1); Q(2, 2) - Q(1, 2), Q(4, 2) - Q(1, 2)];
    temp = Jac_phi' * Jac_phi;
    
    a = (1 / (temp(1, 1) * temp(2, 2) - temp(1, 2) ^ 2)) * temp(2, 2);
    b = (-1 / (temp(1, 1) * temp(2, 2) - temp(1, 2) ^ 2)) * temp(1, 2);
    c = (1 / (temp(1, 1) * temp(2, 2) - temp(1, 2) ^ 2)) * temp(1, 1);
        
    M = (2 * (a + c) + 3 * b) * eye(4);
    M(1:3, 4) = [-2 * c + a; -(a + 3 * b + c); c - 2 * a];
    M(1:2, 3) = [-(a + 3 * b + c); -2 * c + a];
    M(1, 2) = -2 * c + a;

    M(2:4, 1) = [-2 * c + a; -(a + 3 * b + c); c - 2 * a];
    M(3:4, 2) = [c - 2 * a, -(a + 3 * b + c)];
    M(4, 3) = c - 2 * a;
end
