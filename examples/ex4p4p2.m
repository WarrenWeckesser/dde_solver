function [y,ylag] = ex4p4p2
temp = load('ex4p4p2.dat','-ascii');
y = temp(:,1);
ylag = temp(:,2);

plot(y,ylag)
title('ex4p4p2')
xlabel('y(t)');
ylabel('y(t - 2)');
