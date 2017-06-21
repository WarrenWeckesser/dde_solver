function [t,y] = ex4p4p1
temp = load('ex4p4p1.dat','-ascii');
t = temp(:,1);
y = temp(:,2:end);

plot(t,y)
title('ex4p4p1')
legend('S(t)','E(t)','I(t)','R(t)')