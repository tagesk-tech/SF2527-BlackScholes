function [x, u, tau, u_history] = forward_euler(K, r, sigma, T, h, Smax, Nt, plot_limit)
x = 0:h:Smax;
u = max(x - K, 0);

dt = T/Nt;
m = length(x);
plot_m = find(x <= plot_limit, 1, 'last');
tau_all = (0:Nt)*dt;
u_history = zeros(Nt + 1, plot_m);
u_history(1, :) = u(1:plot_m);

for n = 1:Nt
    mnew = m - 1;
    i = 2:mnew;
    diffusion = sigma^2*x(i).^2/(2*h^2);
    advection = r*x(i)/(2*h);

    unew = zeros(1, mnew);
    unew(i) = u(i) + dt*( ...
        diffusion.*(u(i+1) - 2*u(i) + u(i-1)) ...
        + advection.*(u(i+1) - u(i-1)) - r*u(i));

    u(1:mnew) = unew;
    u_history(n + 1, :) = u(1:plot_m);
    m = mnew;
end

x = x(1:m);
u = u(1:m);
tau = tau_all(1:n+1);
u_history = u_history(1:n+1, :);
end
