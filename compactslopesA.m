function A = compactslopesA(h)

n = length( h ) +1;
hmean = mean(h);
A = spalloc( n, n, 3*n-2 ); % sparse( n, n );

% Special formula at the left end
A(1,1) = h(2)*(h(3)+h(2))/(h(2)+h(1))/(h(3)+h(2)+h(1));
A(1,2) = 1.0;

% Symmetric compact finite difference formula on interior nodes
for k=2:n-1
    A(k,k-1) = 4*h(k)^2/(h(k-1)+h(k))^2;
    A(k,k)   = 4;
    A(k,k+1) =  4*h(k-1)^2/(h(k-1)+h(k))^2;
end


% Special formula at the right end
A(n,n-1) = 1.0;
A(n,n)   = (h(n-2))*(h(n-2)+h(n-3))/(h(n-1)+h(n-2))/(h(n-3)+h(n-2)+h(n-1) );

