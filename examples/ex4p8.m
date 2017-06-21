function [t,y] = ex4p8

temp = load('export.dat','-ascii');
t = temp(:,1);
y = temp(:,2:end);

% Scale so as to reproduce Fig. 15.8 of Hairer, Norsett,
% and Wanner book.
yplot = [1e4*y(:,1) 0.5*y(:,2) y(:,3) 10*y(:,4)];
plot(t,yplot)
axis([0 60 -1 15.5])
legend('10^4 V','C/2','F','10m',0)
