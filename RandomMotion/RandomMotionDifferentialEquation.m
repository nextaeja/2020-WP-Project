function dydt = RandomMotionDifferentialEquation(t, y, xi, xi_t, gamma, m)

dydt = zeros(2,1);
xit = interp1(xi_t, xi, t);
dydt(1) = y(2);
dydt(2) = -gamma*y(2) + 500*gamma*xit;
%fprintf('%eps\n', t/1e-12);