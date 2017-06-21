function [t,y] = c2g

temp = load('c2g.dat','-ascii');
t = temp(:,1);
y = temp(:,2:end);

plot(t,y)
title('C2')

