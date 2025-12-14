% Fourier series sketch for:
% f(x) = x,   -L < x < 0
%       2x,   0 < x < L
clear; clc;

L = 1; 
Ns = [3 50 100];
xx = linspace(-L, L, 4000);

% original piecewise f on (-L,L) (one period)
ftrue = (xx<0).*xx + (xx>=0).*2.*xx;

% Fourier coefficients (period 2L):
% f(x) ~ a0/2 + sum_{n>=1} [ a_n cos(n*pi*x/L) + b_n sin(n*pi*x/L) ]
% a0 = L/2
% a_n = L*((-1)^n - 1)/(n^2*pi^2)   (=> a_n = 0 for even n; = -2L/((2k-1)^2*pi^2) for odd n)
% b_n = -3*L*(-1)^n/(n*pi)

a0 = L/2;

figure; hold on;
plot(xx, ftrue, 'k', 'LineWidth', 1.8);  % original f

for N = Ns
    S = (a0/2)*ones(size(xx));
    for n = 1:N
        an = L*((-1)^n - 1)/(n^2*pi^2);
        bn = -3*L*(-1)^n/(n*pi);
        S = S + an*cos(n*pi*xx/L) + bn*sin(n*pi*xx/L);
    end
    plot(xx, S, 'LineWidth', 1.2, 'DisplayName', sprintf('N = %d', N));
end

grid on; box on;
xlabel('x'); ylabel('value');
title('Fourier series partial sums vs. f(x) on [-L,L]');
legend('f(x)', 'N=3','N=50','N=100', 'Location','northwest');
xlim([-L L]);
