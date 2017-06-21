function [tI1,I1,tI2,I2,tI3,I3,tI4,I4] = infectionIoft
temp = load('Ioft.dat');

% Read the number of points in the first integration
% and then the data.  Repeat for the other runs.
e = 0;
npts = temp(e+1,1);
b = e + 2;
e = b + npts-1;
tI1 = temp(b:e,1);
I1 = temp(b:e,2);

npts = temp(e+1,1);
b = e + 2;
e = b + npts-1;
tI2 = temp(b:e,1);
I2 = temp(b:e,2);

npts = temp(e+1,1);
b = e + 2;
e = b + npts-1;
tI3 = temp(b:e,1);
I3 = temp(b:e,2);

npts = temp(e+1,1);
b = e + 2;
e = b + npts-1;
tI4 = temp(b:e,1);
I4 = temp(b:e,2);

figure
plot(tI1,I1,tI2,I2,tI3,I3,tI4,I4)
xlabel('t');
ylabel('I(t)');
title('Hoppensteadt--Jackiewicz Problem')
legend('r(t) = 0.2','r(t) = 0.3','r(t) = 0.4','r(t) = 0.5',0)
axis([0 8 0 10])
