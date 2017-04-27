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
u_minus = 1/2;
u_plus = -1/2;
k = [0:N/2, -N/2+1:-1]*(pi/L);

t = 0;
dt = 1e-2;

u = tanh(-10*x)/2;
v = -5*(sech(10*x)).^2;
v_hat = fft(v);
