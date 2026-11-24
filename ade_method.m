function [x, u, tau, u_history] = ade_method(K, r, sigma, T, h, Smax, Nt, plot_limit)
x = 0:h:Smax;
u = max(x - K, 0);
phi = zeros(size(u));
psi = zeros(size(u));

dt = T/Nt;
m = length(x);
plot_m = find(x <= plot_limit, 1, 'last');
tau = (0:Nt)*dt;
u_history = zeros(Nt + 1, plot_m);
u_history(1, :) = u(1:plot_m);

for n = 1:Nt
    mnew = m - 1;
    phi(1) = 0;

    % Left-to-right sweep, equation (8)
    for i = 2:mnew
        diffusion = sigma^2*x(i)^2/(2*h^2);
        advection = r*x(i)/(2*h);

        numerator = u(i)/dt ...
            + diffusion*(u(i+1) - u(i) + phi(i-1)) ...
            + advection*(u(i+1) - u(i) - phi(i-1)) ...
            - r*u(i)/2;
        denominator = 1/dt + diffusion - advection + r/2;
        phi(i) = numerator/denominator;
    end

    % Right-to-left sweep, equation (9)
    psi(mnew) = phi(mnew);
    for i = mnew-1:-1:2
        diffusion = sigma^2*x(i)^2/(2*h^2);
        advection = r*x(i)/(2*h);

        numerator = u(i)/dt ...
            + diffusion*(psi(i+1) - u(i) + u(i-1)) ...
            + advection*(psi(i+1) + u(i) - u(i-1)) ...
            - r*u(i)/2;
        denominator = 1/dt + diffusion + advection + r/2;
        psi(i) = numerator/denominator;
    end

    psi(1) = 0;
    u(1:mnew) = (phi(1:mnew) + psi(1:mnew))/2; % equation (10)
    u_history(n + 1, :) = u(1:plot_m);
    m = mnew;
end

x = x(1:m);
u = u(1:m);
end
