%% Problem 3.3.1  f(x)=1+3x on 0<x<5
clear; clc; close all
L = 5;  T = 2*L; % period 10
xs = linspace(-15,15,6000);

% odd (sine) extension:
y = mod(xs+L,T) - L; % map to (-L, L]
f_odd  = (y>0).*(1+3*y) + (y<0).*(3*y-1);  % at y=0 or +/-L the series -> 0
f_even = 1 + 3*abs(y); % even (cosine) extension


% Half-range sine: f(x) ~ sum_{n>=1} b_n sin(n*pi*x/L) on (0,L)
% b_n = (2/(n*pi))*(1 - 16*(-1)^n)
bs = @(n) (2./(n*pi)).*(1 - 16*(-1).^n);

% Half-range cosine: f(x) ~ a0/2 + sum a_n cos(n*pi*x/L) on (0,L)
% a0 = (2/L)int_0^L (1+3x) dx = 17
% a_n = 30*((-1)^n - 1)/(n^2*pi^2) -> = -60/(n^2*pi^2) for n odd, 0 for n even
a0 = 17;
ac = @(n) 30*((-1).^n - 1)./( (n.^2)*(pi^2) );

% partial sums and plots
Ns = [5 25 100];

% Sine series
figure; hold on
for N = Ns
    S = zeros(size(xs));
    for n = 1:N
        S = S + bs(n)*sin(n*pi*xs/L);
    end
    plot(xs,S,'LineWidth',1.2,'DisplayName',sprintf('Sine N=%d',N));
end
plot(xs,f_odd,'k','LineWidth',1.5,'DisplayName','odd periodic f');
title('Half-range SINE series of f on [-15,15]'); grid on; legend show

% Cosine series
figure; hold on
for N = Ns
    C = (a0/2)*ones(size(xs));
    for n = 1:N
        C = C + ac(n)*cos(n*pi*xs/L);
    end
    plot(xs,C,'LineWidth',1.2,'DisplayName',sprintf('Cosine N=%d',N));
end
plot(xs,f_even,'k','LineWidth',1.5,'DisplayName','even periodic f');
title('Half-range COSINE series of f on [-15,15]'); grid on; legend show

%% (c) and (d)
pts = [5 9 22 101];

% helper that returns value of odd/even periodic extensions at points,
% with the convention: at y=0 or +/- L the odd series -> 0; even series -> 1+3*|+/- L| = 16
wrap = @(x) mod(x+L,T)-L;
ypts = wrap(pts);

sine_vals   = (ypts==0 | abs(ypts)==L).*0 + ... % endpoints -> 0
              (ypts>0).*(1+3*ypts) + (ypts<0).*(3*ypts-1);

cosine_vals = 1 + 3*abs(ypts); % endpoints -> 16

disp(table(pts.', ypts.', sine_vals.', cosine_vals.', ...
    'VariableNames', {'x','y_in_(-5,5]','SineSeriesValue','CosineSeriesValue'}));
