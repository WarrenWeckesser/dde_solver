function [t,y] = neutvan
temp = load('neutvan.dat','-ascii');
t = temp(:,1);
y = temp(:,2:end);

plot(t,y)
title('neutvan')
