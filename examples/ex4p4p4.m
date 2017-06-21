function [t,y,yp,te,ye,ie] = ex4p4p4
temp = load('ex4p4p4.dat','-ascii');
t = temp(:,1);
y = temp(:,2);
yp = temp(:,3);

plot(y,yp)
title('ex4p4p4--phase plane')
xlabel('y(t)');
ylabel('y''(t)');

temp = load('ex4p4p4extr.dat','-ascii');
ie = temp(:,1);
te = temp(:,2);
ye = temp(:,3);
n1 = find(ie == 1);
x1 = te(n1);
y1 = ye(n1);
n2 = find(ie == 2);
x2 = te(n2);
y2 = ye(n2);

figure
plot(t,y,'k',x1,y1,'rs',x2,y2,'bo')
title('ex4p4p4--time series with extrema marked')
xlabel('Time t');
ylabel('y(t)')