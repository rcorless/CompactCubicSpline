%
% Use the leading principal minors to get the
% positive factors of the totally nonnegative 
% matrix A, using only positive arithmetic
% (assuming the h are all aligned in the same complex
% direction)
%
% Copyright (c) Robert M. Corless October 2022
% MIT License

% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

% TNFactorT 
% Input: a vector or list of mesh widths.  If they all have the same
% complex sign, then the resulting Matrix for compact finite differences
% is totally nonnegative.  
%
% Output: the L and U factors (each bidiagonal) of the compact
% finite difference matrix.
%
% Error conditions: if only three widths are given, the matrix is singular.
% In that case just fit a cubic and differentiate that.  If only four widths
% are given, the matrix is totally nonnegative but the recurrence used here
% fails; in that case, just build the matrix and factor it normally.  This can 
% even be done symbolically, once and for all, although that doesn't seem 
% very useful (the result is complicated).
%
% If widths are given which do *not* have the same sign, the matrix might
% be singular, and this factoring might fail even if the matrix is not 
% singular (just one zero minor would cause it to fail).  In that case 
% use some other method of solution.
%TNFactorT = proc( h::{list,Vector} )
%  local a, k, L, n, q, U;
function [L,U] = TNFactorT( h )
  
  n = length(h);
  if n <=4
    error("This routine needs at least 5 intervals.")
  end
  
  L = sparse(n+1,n+1);
  U = sparse(n+1,n+1);
  a = zeros(1,n+2);
  A = zeros(1,n+1); % ratios a(k)/a(k-1)
  q = zeros(1,n+2);
  Q = zeros(1,n+1); % ratios q(k)/a(k)
  a(1) = h(2)*(h(3)+h(2))/(h(2)+h(1))/(h(3)+h(2)+h(1));
  U(1,1) = a(1);
  U(1,2) = 1;
  L(1,1) = 1;

  a(2) = h(1)*h(2)*h(3)/(h(1)+h(2)+h(3));
  A(2) = a(2)/a(1);
  U(2,2) = 4*h(3)*h(1)/(h(2)+h(1))/(h(3)+h(2));
  U(2,3) = 4*h(1)^2/(h(1)+h(2))^2;
  L(2,1) = 4*(h(2))^2/(h(1)+h(2))^2/U(1,1);
  L(2,2) = 1;
   
  q(2) = h(1)*h(2)^2*h(3)^2*(h(1)+2*h(2)+h(3))/(h(1)+h(2))/(h(1)+h(2)+h(3));
  Q(2) = q(2)/a(2);
  a(3) = h(1)*h(2)^2*h(3)*(h(2)+h(3))/(h(1)+h(2));
  A(3) = a(3)/a(2);
  U(3,3) = 4*A(3)/(h(2)+h(3))^2;
  U(3,4) = 4*h(2)^2/(h(2)+h(3))^2;
  L(3,2) = 4*(h(3))^2/(h(2)+h(3))^2/U(2,2);
  L(3,3) = 1;
    

  for k=4:n
    %q(k-1) = 2*h(k)*h(k-1)*a(k-1) + h(k)^2*q(k-2);
    %a(k) = h(k-1)^2*a(k-1) + q(k-1);
    Q(k-1) = 2*h(k)*h(k-1) + h(k)^2*Q(k-2)/A(k-1);
    A(k) = h(k-1)^2 + Q(k-1);
    U(k,k) = 4*A(k)/(h(k-1)+h(k))^2;
    U(k,k+1) = 4*h(k-1)^2/(h(k-1)+h(k))^2;
    L(k,k-1) = 4*h(k)^2/(h(k-1)+h(k))^2/U(k-1,k-1);
    L(k,k) = 1;
  end
  
  L(n+1,n) = 1/U(n,n);

  % This formula would have to be modified if n=5
  A(n+1) = h(n-1)^2*h(n-2)^2*h(n)*(h(n-2)+2*h(n-1)+h(n))/ ...
            ((h(n-1)+h(n))*(h(n-2)+h(n-1)+h(n))*A(n)*A(n-1)) + ...
            h(n)*h(n-2)*h(n-1)^3*Q(n-3)/ ...
            (A(n)*A(n-1)*A(n-2)* ...
            (h(n-2)+h(n-1)+h(n)));
  U(n+1,n+1) = A(n+1);
  L(n+1,n+1) = 1;
end