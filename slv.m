%
% Direct solver (bypasses singularity tests)
%
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


function d = slv(L,U,y)
  n = length(y)-1;
  tmp = zeros(n+1,1);
  tmp(1) = y(1);
  for k = 2:n+1
    tmp(k) = -L(k,k - 1) * tmp(k - 1) + y(k);
  end
  d = zeros(n+1,1);
  d(n + 1) = tmp(n + 1) / U(n + 1,n + 1);
  for k = n:-1:1
    d(k) = (-U(k,k + 1) * d(k + 1) + tmp(k)) / U(k,k);
  end
 