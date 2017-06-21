function [t,y] = volinteq

temp = load('volinteq.dat','-ascii');
t = temp(:,1);
y = temp(:,2:end);

plot(t,y)
title('Volterra Integral Equation')
