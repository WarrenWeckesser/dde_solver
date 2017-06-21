function [t,y] = b2g

temp = load('b2g.dat','-ascii');
t = temp(:,1);
y = temp(:,2:end);

plot(t,y)
title('B2')

