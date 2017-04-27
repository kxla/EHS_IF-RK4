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

t_0 = 0;
dt = 1e-2;
t_max = 10;

u = tanh(-10*x)/2;
v = -5*(sech(10*x)).^2;
v_hat = fft(v);
V_hat = exp(mu*1i.*(k.^3).*t_0);

for t = t_0:dt:t_max
  a = dt .* (-1i*c_0.*k .* exp(1i*mu.*(k.^3).*t) .* fft((u.^p).*v)-nu.*(k.^2).*V_hat);
  b = dt .* (-1i*c_0.*k .* exp(1i*mu.*(k.^3).*(t+dt/2)) .* fft((u.^p).*v)-nu.*(k.^2).*(V_hat+a./2));
  c = dt .* (-1i*c_0.*k .* exp(1i*mu.*(k.^3).*(t+dt/2)) .* fft((u.^p).*v)-nu.*(k.^2).*(V_hat+b./2));
  d = dt .* (-1i*c_0.*k .* exp(1i*mu.*(k.^3).*(t+dt)) .* fft((u.^p).*v)-nu.*(k.^2).*(V_hat+c));
  V_hat = V_hat + (1/6).*(a + 2.*b + 2.*c + d);
end
