function [t,y,te,ye,ie] = ex4p7

temp = load('export.dat','-ascii');
t = temp(:,1);
y = temp(:,2:end);

temp = load('events.dat','-ascii');
ie = temp(:,1);
te = temp(:,2);
ye = temp(:,3:end);

plot(t,y,te,ye,'ro')


