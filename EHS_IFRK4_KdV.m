% EHS IF-RK4 KdV
% Method given by El, Hoefer, Shearer
% Solves the Riemann problem u_t + c_0 u^p u_x = nu u_xx + mu u_xxx

tic

nplots = 50;

c_0 = 1;
p = 1;
nu = 1;
mu = 1;
L = 400;
N = 4000;
dx = 2*L/N;
x = (dx-L):dx:L;
u_m = 1/2;
u_p = -1/2;
k = [0:N/2, -N/2+1:-1]*(pi/L);

t = 0;
dt = 1e-3;
t_max = 10;

u = tanh(-10*x)/2+(tanh(10*(x+1))-1)/4+(tanh(10*(x-1))+1)/4;
v = -5*(sech(10*x)).^2;
V_hat = fft(v);

tdata = zeros(nplots+1,1);
uu = u;zeros(N,t_max/dt);

for k3=1:nplots

  for index=1:round(t_max/dt/nplots)
  
    a = dt .* (-1i*c_0.*k .* exp(1i*mu.*(k.^3).*t) .* fft((u.^p).*v)-nu.*(k.^2).*V_hat);
    b = dt .* (-1i*c_0.*k .* exp(1i*mu.*(k.^3).*(t+dt/2)) .* fft((u.^p).*v)-nu.*(k.^2).*(V_hat+a./2));
    c = dt .* (-1i*c_0.*k .* exp(1i*mu.*(k.^3).*(t+dt/2)) .* fft((u.^p).*v)-nu.*(k.^2).*(V_hat+b./2));
    d = dt .* (-1i*c_0.*k .* exp(1i*mu.*(k.^3).*(t+dt)) .* fft((u.^p).*v)-nu.*(k.^2).*(V_hat+c));
  
    V_hat = V_hat + (1/6).*(a + 2.*b + 2.*c + d);
    v_hat = exp(-1i*mu.*(k.^3).*t).*V_hat;
    v = ifft(exp(-1i*mu*(k.^3)*t).*V_hat);
  
    for j = 2:N
      u(j) = sum(v_hat(2:end)./(1i*k(2:end)).*exp(1i*k(2:end)*x(j)));
    end
  
    u = u - trapz(x,x.*v)/(2*L) + (x+L)*(u_p-u_m)/(2*L)+u_m;
        
    t = t + dt;
  
  end
  
  tdata(k3+1,1) = t;
  uu(k3+1,:) = u;
  
end

figure
waterfall(x,tdata,real(uu)), view(0,70),
xlim([-L,L]);
ylim([0,t_max]);
grid off

toc
