function [t,y] = b2

temp = load('export.dat','-ascii');
t = temp(:,1);
y = temp(:,2:end);

plot(t,y)
title('B2')

