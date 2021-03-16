function M = raideur_quadrangle(Q)    
    Jac_phi = [Q(2, 1) - Q(1, 1), Q(4, 1) - Q(1, 1); Q(2, 2) - Q(1, 2), Q(4, 2) - Q(1, 2)];
    temp = Jac_phi' * Jac_phi;
    
    a = (1 / (temp(1, 1) * temp(2, 2) - temp(1, 2) ^ 2)) * temp(2, 2);
    b = (-1 / (temp(1, 1) * temp(2, 2) - temp(1, 2) ^ 2)) * temp(1, 2);
    c = (1 / (temp(1, 1) * temp(2, 2) - temp(1, 2) ^ 2)) * temp(1, 1);
    
    M = [
        2*a+3*b+2*c, -2*a+c, -a-3*b-c, a-2*c;
        -2*a+c, 2*a-3*b+2*c, a-2*c, -a+3*b-c;
        -a-3*b-c, a-2*c, 2*a+3*b+2*c, -2*a*c;
        a-2*c, -a+3*b-c, -2*a+c, 2*a-3 * b+2*c; 
    ];
    
    M = det(Jac_phi) .* M;
end
