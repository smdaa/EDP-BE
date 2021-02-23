function M=raideur_triangle(A)
    % A= [A1(1) , A1(2) ; A2(1) , A2(2) ; A3(1), A3(2) ]
    alpha=det([ (A(2,1)-A(1,1)) , (A(3,1)-A(1,1)) ; (A(2,2)-A(1,2)) , (A(3,2)-A(1,2)) ]);

    M=zeros(3,3);
    
    for i=1:3 
        for j=1:3
            grad_nj=(1/alpha)*[ (A(mod(j+1,3),2)-A(mod(j+2,3),2)) ; (A(mod(j+2,3),1)-A(mod(j+1,3),1))];
            grad_ni=(1/alpha)*[ (A(mod(i+1,3),2)-A(mod(i+2,3),2)) ; (A(mod(i+2,3),1)-A(mod(i+1,3),1))];
            M(i,j)=(alpha/2)*grad_ni'*grad_nj;
        end
    end
    
end