function [t,y] = exsd1

temp = load('exsd1.dat','-ascii');
t = temp(:,1);
y = temp(:,2:end);

plot(t,y)
title('exsd1')
