function ex06ch4
global c D G n
c = 0.5; D = 5; n = 100;
% x = y(1), y = y(2), lambda = y(3), 
% I_1 = y(4), I_2 = y(5), I_3 = y(6).
y0 = [0.8*n; 0.2*n; 0; 0; 0; 0];
for G = [0.1 1 2]
    fprintf('\nFor G = %4.1f, \n',G)
    lambda_s = c*D - G;
    x_s = G*n/(lambda_s + G);
    y_s = n - x_s;
    fprintf('Analytical x_s = %4.1f, y_s = %4.1f, lambda_s = %4.1f.\n',...
            x_s,y_s,lambda_s);
    sol = dde23(@odes,[],y0,[0, D]);
    sol = dde23(@ddes,D,sol,[D, 4*D]);
    fprintf('Computed   x_s = %4.1f, y_s = %4.1f, lambda_s = %4.1f.\n',...
             sol.y(1:3,end));
end
%============================================================
function dydt = odes(t,y,Z)
global c D G n
dydt = zeros(6,1);
dydt(1) = - y(1)*y(3) + G*y(2);
dydt(2) = - dydt(1);
dydt(4) = exp(G*t)*y(2);
dydt(5) = t*exp(G*t)*y(1)*y(3);
dydt(6) = exp(G*t)*y(1)*y(3);
dydt(3) = (c/n)*exp(-G*t)*((dydt(4)+dydt(5))-G*(y(4)+y(5)));

function dydt = ddes(t,y,Z)
global c D G n
dydt = zeros(6,1);
dydt(1) = - y(1)*y(3) + G*y(2);
dydt(2) = - dydt(1);
dydt(4) = exp(G*t)*y(2) - exp(G*(t - D))*Z(2);
dydt(5) = D*exp(G*t)*y(1)*y(3) - y(6);
dydt(6) = exp(G*t)*y(1)*y(3) - exp(G*(t - D))*Z(1)*Z(3);
dydt(3) = (c/n)*exp(-G*t)*((dydt(4)+dydt(5))-G*(y(4)+y(5)));
