%% Test condition numbers of compact finite difference matrices
clear all
ns = Fibonacci(15);
ii = find(ns>100);
ns = ns(ii);
cnds = zeros(1,length(ns));
fn = @(t) 1.0./(1+25*t.^2);
dfn = @(t) -50*t./(1+25*t.^2).^2;
for k=1:length(ns)
    n = ns(k);
    tau = cos(pi*(n:-1:0)/n);
    A = compactslopesA(diff(tau));
    cnds(k) = cond(full(A));
end
figure(1),semilogx( ns, cnds, 'ko', ns, cnds(end)*(ns/ns(end)).^0, 'k-' )
set(gca,'fontsize',16)
xlabel( 'Dimension','fontsize',16 ), ylabel( 'Condition number Cheby','fontsize',16 )
%
for k=1:length(ns)
    n = ns(k);
    tau = linspace(-1,1,n+1);
    A = compactslopesA(diff(tau));
    cnds(k) = cond(full(A));
end
figure(2),semilogx( ns, cnds, 'ko', ns, cnds(end)*(ns/ns(end)).^0, 'k-' )
set(gca,'fontsize',16)
xlabel( 'Dimension','fontsize',16 ), ylabel( 'Condition number equal','fontsize',16 )
%
%
%
% mu = 110;
% [t,y] = vdpode(mu);
% f = y(:,1)';
% df = vcompact4(t,f);
% yy = [f;df'];
% u = RefineMesh(t',25);
% [v,vp,vpp] = cubichermite( t', yy, u );
% df = y(:,2)';
% h = diff(t);
% figure,semilogy(t(2:end), h, 'k.' )
% ddf = vcompact4(t,df);
% yy = [df;ddf'];
% [vp,vpp,vppp] = cubichermite( t', yy, u );
% res = vpp - mu*(1-v.^2).*vp + v;
% plot( u, res/max(abs(vpp)), 'k' )
