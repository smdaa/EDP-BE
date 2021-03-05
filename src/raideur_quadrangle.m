function [M, alpha] = raideur_quadrangle(C)
    alpha=det([ C(2,1) - C(1,1) C(3,1) - C(1,1) ;C(2,2) - C(1,2) C(3,2) - C(1,2)]);
    Jac_phi = [C(2, 1) - C(1, 1) C(4, 1) - C(1, 1); C(2, 2) - C(1, 2) C(4, 2) - C(1, 2)];
    temp = pinv(Jac_phi' * Jac_phi);
    a = temp(1, 1);
    b = temp(1, 2);
    c = temp(2, 1);
        
    M = (2 * (a + c) + 3 * b) * eye(4);
    M(1:3, 4) = [-2 * c + a; -(a + 3 * b + c); c - 2 * a];
    M(1:2, 3) = [-(a + 3 * b + c); -2 * c + a];
    M(1, 2) = -2 * c + a;
    
    M = (M + M.') / 2;
end