function [t,y] = secdelay

temp = load('secdelay.dat','-ascii');
t = temp(:,1);
y = temp(:,2:end);

plot(t,y)
title('secdelay')
