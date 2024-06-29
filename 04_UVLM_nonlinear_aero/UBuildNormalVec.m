function n = UBuildNormalVec(M,N,x,y,z)

grid = cell(M+1,N+1);
for i=1:M+1
    for j=1:N+1
        grid{i,j} = [x(i,j) y(i,j) z(i,j)];
    end
end

n  = zeros(3,M*N);
for i = 1:M
    for j = 1:N
        k = sub2ind([M N],i,j);
        A = grid{i+1,j+1} - grid{i,j};
        B = grid{i,j+1} - grid{i+1,j};
        n(:,k) = cross(A,B)/norm(cross(A,B));
    end
end
%n = n(:);
end

