%
% Test of TNFactorT and of applyB
%
f = @(x) 1.0./(1+25*x.^2); % Runge function
df = @(x) -50.0*x./(1+25*x.^2).^2; % its analytic derivative (by hand)

tau = linspace(-1,1,101);
fv = f(tau);
h = diff(tau); % all positive, all equal

[L,U] = TNFactorT( h );
y = applyB( h, fv );

d = slv( L, U, y );

refd = df(tau); % Row vector, whereas d is column vector

error = d - refd';
figure(1), plot( tau, error, 'k.' )
grid on

% Now do a varying grid
N = 1025;
tau = cos(pi*(N:-1:0)/N);  % Chebyshev points
fv = f(tau);
h = diff(tau); % all positive, all equal

[L,U] = TNFactorT( h );
y = applyB( h, fv );

d = slv( L, U, y );
refd = df(tau); % Row vector, whereas d is column vector

err = d - refd';
figure(2), plot( tau, abs(err), 'k.' )
grid on

