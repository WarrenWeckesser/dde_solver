function [t1,y1,t2,y2,t3,y3,t4,y4] = infection8
temp = load('infection.dat');
neqn = temp(1);

% Read the number of points in the first integration
% and then the data.  Repeat for the other runs.
e = 1;
npts = temp(e+1);
b = e + 2;
e = b + npts*(neqn+1)-1;
M = reshape(temp(b:e),neqn+1,npts)';
t1 = M(:,1);
y1 = M(:,2:neqn+1);

npts = temp(e+1,1);
b = e + 2;
e = b + npts*(neqn+1)-1;
M = reshape(temp(b:e),neqn+1,npts)';
t2 = M(:,1);
y2 = M(:,2:neqn+1);

npts = temp(e+1,1);
b = e + 2;
e = b + npts*(neqn+1)-1;
M = reshape(temp(b:e),neqn+1,npts)';
t3= M(:,1);
y3 = M(:,2:neqn+1);

npts = temp(e+1,1);
b = e + 2;
e = b + npts*(neqn+1)-1;
M = reshape(temp(b:e),neqn+1,npts)';
t4 = M(:,1);
y4 = M(:,2:neqn+1);

figure
plot(t1,y1(:,2),t2,y2(:,2),t3,y3(:,2),t4,y4(:,2))
xlabel('t');
ylabel('S(t)');
title('Hoppensteadt--Jackiewicz Problem')
legend('r(t) = 0.2','r(t) = 0.3','r(t) = 0.4','r(t) = 0.5',0)
