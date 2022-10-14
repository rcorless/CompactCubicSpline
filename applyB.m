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
%
% To compute compact finite differences from a vector of values 
% of a function f, one has to solve the linear system Td = Bf
% where B is a particular sparse matrix.
%
% The routine "applyB" builds the elements of B from the vector
% of widths and applies it to the vector of function values.
% After that, the L and U factors from the previous routine can
% be used to find the derivatives by solving
% 
% Lv = y (= B*f)
% Ud = v
%
% in succession.  This is done in the worksheet of this Maple workbok.
%
function y = applyB( h, f )
  n = length(h);
  y = zeros(n+1,1);
  a = zeros(1,4);
  a(1) = -(h(2)+h(3))*h(2)/h(1)*(4*h(1)^2+6*h(1)*h(2)+3*h(1)*h(3)+...
           2*h(2)^2+2*h(2)*h(3))/(h(1)+h(2))^2/(h(1)+h(2)+h(3))^2;
  a(2) = -1/h(1)*(2*h(1)*h(2)+h(1)*h(3)-2*h(2)^2-2*h(2)*h(3))/h(2)/(h(2)+h(3));
  a(3) = (h(2)+h(3))/h(2)*h(1)^2/(h(1)+h(2))^2/h(3);
  a(4) = -h(1)^2*h(2)/h(3)/(h(2)+h(3))/(h(1)+h(2)+h(3))^2;
  
  y(1) = a(1)*f(1) + a(2)*f(2) + a(3)*f(3) + a(4)*f(4);
  
  for k = 2:n 
    y(k) = -8*h(k)^2*(2*h(k-1)+h(k))*f(k-1)/(h(k-1)*(h(k-1)+h(k))^3) ...
            -8*(h(k-1)-h(k))*f(k)/(h(k-1)*h(k)) ...
            +8*h(k-1)^2*(h(k-1)+2*h(k))*f(k+1)/(h(k)*(h(k-1)+h(k))^3);
  end
  
  a(1) = h(n-1)*(h(n-2)+h(n-1))*(4*h(n)^2+3*h(n)*h(n-2)+6*h(n)*h(n-1) ...
               +2*h(n-2)*h(n-1)+2*h(n-1)^2)/(h(n)+h(n-1)+h(n-2))^2/ ...
               (h(n)+h(n-1))^2/h(n);
  a(2) = (h(n)*h(n-2)+2*h(n)*h(n-1)-2*h(n-2)*h(n-1)-2*h(n-1)^2)/h(n) ...
           /h(n-1)/(h(n-2)+h(n-1));
  a(3) = -(h(n-2)+h(n-1))*h(n)^2/h(n-2)/h(n-1)/(h(n)+h(n-1))^2;
  a(4) = h(n-1)*h(n)^2/h(n-2)/(h(n-2)+h(n-1))/(h(n)+h(n-1)+h(n-2))^2;
  
  y(n+1) = a(1)*f(n+1) + a(2)*f(n) + a(3)*f(n-1) + a(4)*f(n-2);

end 
