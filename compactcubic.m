function [v,vp,vpp] = compactcubic(xin,yin,uin)
%COMPACTCUBIC  compact cubic spline alternative.
%  [v,vp,vpp] = compactcubic(x,y,u) finds the piecewise cubic interpolatory
%  function S(x), with S(x(j)) = y(j), and returns v(k) = S(u(k)) and its
%  derivatives.  Adapted by RMC from splinetx.m by Cleve Moler from 
%  Numerical Computing with Matlab
%
%  See also vcompact4.m

   n = length(xin);
   x = zeros(1,n);
   x(:) = xin;
   y = zeros(1,n);
   y(:) = yin;
   u = zeros(1,length(uin));
   u(:) = uin;

   h = diff(x);
   if any(h<=0) 
       error("nodes must be strictly increasing" )
   end
   d = compactslopes(h,y).';
   delta = diff(y)./h;
   ms = (d(1:end-1)+d(2:end)); % 2Mean slopes

%  Piecewise polynomial coefficients

   c = (3*delta - ms - d(1:end-1))./h;
   b = (ms - 2*delta )./h.^2; % Basically the same formulas

%  Find subinterval indices k so that x(k) <= u < x(k+1)

   k = ones(size(u));
   for j = 2:n-1
      k(x(j) <= u) = j;
   end

%  Evaluate spline
   % and its derivatives 
   s = u - x(k);
   v = y(k) + s.*(d(k) + s.*(c(k) + s.*b(k)));
   vp = d(k) + s.*(2*c(k) + s.*(3*b(k)));
   vpp= 2*c(k) + s.*(6*b(k));


% -------------------------------------------------------

function df = compactslopes(h,f)

n = length( f );

A = sparse( n, n );
B = sparse( n, n );

% Special formula at the left end
A(1,1) = h(2)*(h(3)+h(2))/(h(2)+h(1))/(h(3)+h(2)+h(1));
A(1,2) = 1.0;
% Unless carefully written, these are quite susceptible to rounding errors.
B(1,1) = (4*h(1)^2+6*h(1)*h(2)+3*h(1)*h(3)+2*h(2)^2+2*h(2)*h(3)) ...
         *h(2)*(h(3)+h(2))/(h(1)+h(2))^2/(h(1)+h(2)+h(3))^2/h(1);
B(1,2) = 1/h(1)*((-2*h(2)+h(1))*h(3)+2*h(2)*(-h(2)+h(1))) ...
                 /h(2)/(h(2)+h(3));
B(1,3) = -h(1)^2*(h(3)+h(2))/(h(2)+h(1))^2/(h(2))/h(3);
B(1,4) =  h(1)^2*h(2)/(h(3)+h(2)+h(1))^2/(h(3)+h(2))/h(3);

% Symmetric compact finite difference formula on interior nodes
for k=2:n-1
    A(k,k-1) = 4*h(k)^2/(h(k-1)+h(k))^2;
    A(k,k)   = 4;
    A(k,k+1) =  4*h(k-1)^2/(h(k-1)+h(k))^2;
    B(k,k-1) = 8*h(k)^2*(2*h(k-1)+h(k))/h(k-1)/(h(k-1)+h(k))^3;
    B(k,k)   = -8*(-h(k-1)+h(k))/h(k-1)/h(k);
    B(k,k+1) = -8*h(k-1)^2*(h(k-1)+2*h(k))/(h(k-1)+h(k))^3/h(k);
end


% Special formula at the right end
A(n,n-1) = 1;
A(n,n)   = (h(n-2))*(h(n-2)+h(n-3))/(h(n-1)+h(n-2))/(h(n-3)+h(n-2)+h(n-1) );
% As above, unless care is used, these are quite susceptible to rounding errors.
B(n,n-3) = -h(n-2)*h(n-1)^2/(h(n-3)+h(n-2)+h(n-1))^2/(h(n-3)+h(n-2))/h(n-3);
B(n,n-2) = h(n-1)^2/(h(n-2)+h(n-1))^2/h(n-2)/h(n-3)*(h(n-2)+h(n-3));
B(n,n-1) = 1/h(n-1)*((2*h(n-2)-h(n-1))*h(n-3)+2*h(n-2)*(h(n-2)-h(n-1))) ...
           /h(n-2)/(h(n-3)+h(n-2));
B(n,n) = -(4*h(n-1)^2+(6*h(n-2)+3*h(n-3))*h(n-1)+2*h(n-2)*(h(n-3)+h(n-2))) ...
           /h(n-1)/(h(n-2)+h(n-1))^2/(h(n-3)+h(n-2)+h(n-1))^2 ...
           *(h(n-2))*(h(n-2)+h(n-3));


% Solve the tridiagonal system for the unknown derivatives
df = -A\(B*f(:));
