%% Test the order of accuracy of the compact cubic spline on uniform and Chebyshev grids
%
% Spoiler: uniform beats Chebyshev (but Chebyshev does ok)
% RMC Jan 2019

ns = Fibonacci(22);
ns = ns(5:end);
errs = zeros(1,length(ns));
derrs = zeros(1,length(ns));
errsC = zeros(1,length(ns));
derrsC = zeros(1,length(ns));
Runge = @(x) 1.0./(1+25*x.^2);
dRunge = @(t) -50*t./(1+25*t.^2).^2;
u = linspace(-1,1,2019);
Ru  =  Runge(u);
dRu = dRunge(u);
for k=1:length(ns)
    n = ns(k);
    tau = linspace(-1,1,n+1);
    f = Runge(tau);
    [v,vp,vpp] = compactcubic( tau, f, u );
    errs(k) = norm( v - Ru, inf );
    derrs(k) = norm( vp - dRu, inf );
end

for k=1:length(ns)
    n = ns(k);
    tau = cos(pi-pi*(0:n)/n);
    f = Runge(tau);
   [v,vp,vpp] = compactcubic( tau, f, u );
    errsC(k) = norm( v - Ru, inf );
    derrsC(k) = norm( vp - dRu, inf );
end
close all
figure(1)
equif = loglog( ns, errs, 'ko' );
xlabel('n')
ylabel('errors (equispaced grid)')
hold on
equifit = loglog( ns, ns.^(-4)*errs(end)*ns(end)^4, 'k-' );
equid = loglog( ns, derrs, 'kx');
equidfit = loglog( ns, ns.^(-3)*derrs(end)*ns(end)^3, 'k-.'  );
hold off
legend([equif equifit equid equidfit],{'function', 'O(1/n^4)', 'derivative', 'O(1/n^3)' }, 'NumColumns', 2)

figure(2)
chebyf = loglog( ns, errsC, 'ko');
xlabel('n')
ylabel('errors (Chebyshev grid)')
hold on
chebyfit = loglog( ns, ns.^(-4)*errsC(end)*ns(end)^4, 'k-'  );
chebyd = loglog( ns, derrsC, 'kx');
chebydfit = loglog( ns, ns.^(-3)*derrsC(end-3)*ns(end-3)^3, 'k-.' );
legend([chebyf chebyfit chebyd chebydfit], {'function', 'O(1/n^4)', 'derivative', 'O(1/n^3)'},  'NumColumns', 2)