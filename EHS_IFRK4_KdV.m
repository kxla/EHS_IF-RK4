% EHS IF-RK4 KdV
% Integrating Factor Method
% Solves u_t + c_0 u u_x = -c_1 u_x + mu u_xxx

tic;

nplots = 50;

c_0 = 1;
c_1 = 2;
mu = 1;
N = 2^10;
Left = 50;
k = [0:N/2, -N/2+1:-1]*(pi/Left);
L = -1i*k*c_1-1i*(k.^3)*mu;
dx = 2*Left/N;
x = (dx-Left):dx:Left;
alpha = 1;
xi = 1;

t = 0;
dt = 5e-3;
t_max = 10;

%{
syms y;
func(y) = tanh(-5*y)/2;
df = diff(func,y);
u = double(subs(func,y,x));
u_m = u(1);
u_p = u(end);
v = double(subs(df,y,x));
%}

u = alpha*sech(xi*x).^2;
u_hat = fft(u);
v = exp(-L*t).*u_hat;

tdata = zeros(nplots+1,1);
uu = u;zeros(N,t_max/dt);

filter = floor(N/3);
if mod(filter,2)==0
    filter = filter - 1;
end

for k3=1:nplots

  for index=1:round(t_max/dt/nplots)
      
    % filter to combat aliasing
    %V_hat(N/2 + 1 - (filter-1)/2:N/2 + 1 + (filter-1)/2) = 0;
    
    a = dt * exp(-L*t).*((-c_0/2)*1i*k).*fft((ifft(exp(L*t).*v)).^2);
    b = dt * exp(-L*(t+(dt/2))).*((-c_0/2)*1i*k).*fft((ifft(exp(L*(t+(dt/2))).*(v+(a/2)))).^2);
    c = dt * exp(-L*(t+(dt/2))).*((-c_0/2)*1i*k).*fft((ifft(exp(L*(t+(dt/2))).*(v+(b/2)))).^2);
    
    t = t+dt;
    
    d = dt * exp(-L*t).*((-c_0/2)*1i*k).*fft((ifft(exp(L*t.*(v+c)))).^2);
    
    v = v + (1/6).*(a+2.*b+2.*c+d);
    
  end
  
  tdata(k3+1,1) = t;
  uu(k3+1,:) = ifft(exp(L*t).*v);
  
end

figure
waterfall(x,tdata,real(uu)), view(0,70),
xlim([-Left,Left]);
ylim([0,t_max]);
grid off

toc
