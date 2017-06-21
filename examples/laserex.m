function [t,y] = laserex
% Import the output of laserex.f90 and plot it as done 
% in T.W. Carr's solution with DDE23. 

temp = load('laserex.dat','-ascii');
t = temp(:,1);
y = temp(:,2:end);

delay = 20;

% The lagged solution y(t - delay) is to be plotted against
% y(t).  It is assumed that for a positive integer 'refine',
% the solution has been computed at the points
%  
%   spacing = delay/refine;
%   nout = 1 + round( (tfinal - t0)/spacing );
%   t = zeros(nout,1);
%   for i = 1:nout
%     t(i) = t0 + (i-1)*spacing
%   end

spacing = t(2) - t(1);
refine = delay / spacing;
tplot = t(refine+1:end);
yplot = y(refine+1:end,:);
ylag  = y(1:end-refine,:);
yratio = ylag(:,1) ./ yplot(:,1);

% Plot the field and inversion.
figure(1)
subplot(2,1,1) 
plot(tplot,yplot(:,1:2),'b',tplot,ylag(:,1),'g--');
subplot(2,1,2)
plot(tplot,yplot(:,3));

% Plot the field and inversion.
figure(2)
plot(tplot,yplot(:,1:2));

figure(3)
plot(tplot,yplot(:,1),'b:',tplot,ylag(:,1),'g:',tplot,yratio);
axis([200 400 0 100])
