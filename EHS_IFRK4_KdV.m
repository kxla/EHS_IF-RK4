% EHS IF-RK4 KdV
% Method given by El, Hoefer, Shearer
% Solves the Riemann problem u_t + c_0 u^p u_x = nu u_xx + mu u_xxx

c_0 = 1;
p = 1;
nu = 1;
mu = 1;
L = 1;
N = 32;
dx = 2*L/N;
x = (dx-L):dx:L;
u_m = 1/2;
u_p = -1/2;
