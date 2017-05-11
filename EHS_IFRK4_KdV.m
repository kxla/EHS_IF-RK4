% EHS IF-RK4 KdV
% Method given by El, Hoefer, Shearer
% Solves the Riemann problem u_t + c_0 u u_x = -c_1 u_x + mu u_xxx

tic;

nplots = 50;

c_0 = 1;
c_1 = 2;
% p = 1;
% nu = 1;
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
dt = 5e-5;
t_max = 0.5;

%{
syms y;
func(y) = tanh(-5*y)/2;
df = diff(func,y);
u = double(subs(func,y,x));
u_m = u(1);
u_p = u(end);
v = double(subs(df,y,x));

V_hat = fft(v);
%}

u = alpha*sech(xi*x).^2;
v = exp(-L*t).*u;

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
    
    %uv_hat = fft(u.*v);
  
    % a = dt .* (-1i*c_0.*k .* exp(1i*mu.*(k.^3).*t) .* uv_hat-nu.*(k.^2).*V_hat);
    % b = dt .* (-1i*c_0.*k .* exp(1i*mu.*(k.^3).*(t+dt/2)) .* uv_hat-nu.*(k.^2).*(V_hat+a./2));
    % c = dt .* (-1i*c_0.*k .* exp(1i*mu.*(k.^3).*(t+dt/2)) .* uv_hat-nu.*(k.^2).*(V_hat+b./2));
    
    a = dt * exp(-L*t).*((-c_0/2)*1i*k).*fft((ifft(exp(L*t).*v)).^2);
    b = dt * exp(-L*(t+(dt/2))).*((-c_0/2)*1i*k).*fft((ifft(exp(L*(t+(dt/2))).*(v+(a/2)))).^2);
    c = dt * exp(-L*(t+(dt/2))).*((-c_0/2)*1i*k).*fft((ifft(exp(L*(t+(dt/2))).*(v+(b/2)))).^2);
    d = dt * exp(-L*(t+dt)).*((-c_0/2)*1i*k).*fft((ifft(exp(L*(t+dt).*(v+c)))).^2);
    
    t = t+dt;
    
    % d = dt .* (-1i*c_0.*k .* exp(1i*mu.*(k.^3).*(t)) .* uv_hat-nu.*(k.^2).*(V_hat+c));
  
    % V_hat = V_hat + (1/6).*(a + 2.*b + 2.*c + d);
    % v_hat = exp(-1i*mu.*(k.^3).*t).*V_hat;
    
    % v = ifft(v_hat);
    
    %u0 = u(N/2);
    
    %u = [u0,ifft(v_hat(2:end)./(1i*k(2:end)))] - trapz(x,x.*v)/(2*L) + (x+L)*(u_p-u_m)/(2*L)+u_m;
    
    v = v + (1/6).*(a+2.*b+2.*c+d);
    
    
    %{
    u = (exp(1i*[x(1:N/2-1),x(N/2+1:end)]'*k(2:end))*(v_hat(2:end)./(1i*k(2:end)))')';
    u = [u(1:N/2-1),u0,u(N/2:end)] - trapz(x,x.*v)/(2*L) + (x+L)*(u_p-u_m)/(2*L)+u_m;
    %}
  end
  
  u = ifft(exp(L*t).*v);
  
  tdata(k3+1,1) = t;
  uu(k3+1,:) = u;
  
end

figure
waterfall(x,tdata,real(uu)), view(0,70),
xlim([-Left,Left]);
ylim([0,t_max]);
grid off

toc
