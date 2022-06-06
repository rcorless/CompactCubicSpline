%% Condition numbers of compact A for random meshes
close all
ns = 10946; 
ns = ns(end:end);
m = 10000;
condsR = zeros(m,length(ns));
rng('default')
for k=1:length(ns)
  n = ns(k); % Only one in this demo loop
  for j=1:m
    h = rand(1,n);
    A = compactslopesA(h);
    condsR(j,k) = condest(A);
  end
end

figure(2), histogram( log(condsR(:,1))/log(10),
           'DisplayStyle','stairs', 'EdgeColor', 'Black');
set(gca,'fontsize',24)
ylabel(sprintf('frequency/%d',m),'fontsize',16)
xlabel('log10(condest(A))','fontsize',16)
hold on
plot( log(ns)/log(10)*[1,1], [0,1000], 'k')
title(sprintf('Distribution of condest(A), n=%d',ns))