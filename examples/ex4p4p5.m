function [t,y] = ex4p4p5
temp = load('ex4p4p5.dat','-ascii');
t = temp(:,1);
y = temp(:,2:end);

plot(y(:,1),y(:,2))
title('ex4p4p5')
xlabel('\theta(t)')
ylabel('\theta''(t)')


