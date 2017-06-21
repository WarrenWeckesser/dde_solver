function [t,y] = ex4p4p3
temp = load('ex4p4p3.dat','-ascii');
t = temp(:,1);
y = temp(:,2:end);

% Scale so as to reproduce Fig. 4.3 of the SGT book.
plot(t,[y(:,1:2) 100*y(:,3)])
title('ex4p4p3')
legend('x(t)','y(t)','100\lambda(t)',0)
axis([0 20 0 160])
