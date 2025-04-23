% ASE324L Lab 8
% Anna Ring - aer3965
T = [23, -23.6, 68.4, -49.5, 56, 40.3, 5.4]';
dE = [56.03, 4.92, 89.01, 5.65, 75.16, 65.107, 30.753]';
EVals = [T, dE];
EValsSort = EVals.*0;
[EValsSort(:, 1), idx] = sortrows(EVals(:, 1));
EValsSort(:, 2) = EVals(idx, 2);
ftLbToNM = 1.3558;
clear EVals EValsUnsort T dE idx

EValsSort(:, 2) = ftLbToNM.*EValsSort(:, 2);

figure;
hold on;
plot(EValsSort(:, 1), EValsSort(:, 2), ':o', "LineWidth", 1.5)
xlabel('Temperature (C)')
ylabel('Energy Release (Nm)')
title("Temperature vs Energy Released")
hold off;

%Q3
c0 = 10/100;
cf = 20/100;
H = 1/200;
B = 0.5/100;
G = 30;
E = 100* 10^9;
cArray = linspace(c0, cf, 100);

PArray = sqrt( (G.*E.*H.^3 .* B.^2) ./ (12.*cArray.^2));
DArray = sqrt( (4.*G.*cArray.^4) ./ (3.*E.*H.^3));
figure;
hold on;
plot(DArray, PArray, 'LineWidth',1.5)
xlabel('Displacement (m)')
ylabel('Force (N)')
title('Force vs Displacement Curve for C 10 to 20 cm')
grid on;
hold off;

W = (10/100);
b = 1/100;
P = 100*1000;
ai = 5e-6;
%af = 7.22e-3;
Kic = 11*10^6;

af = (1/pi).*((Kic.*W.*b) / (0.73*P))^2;

fun = @(x) 1./((0.73.*P.*sqrt(pi.*x)./(W.*b*10^6)).^20);
fun2 = @(x) (0.73.*P.*sqrt(pi.*x)./(W.*b.*10^6)).^(-20);

t = integral(fun2, ai, af);
tf = t./2;

M = 200;
a05 = 1.5/1000;
A = 5e-4;
n = 12;
K5 = 6e6;

% A = A/1000;
% A = A*(24*60*60);

af5 = (3.975.*M.*10.^(-6)).^(2/3);

fun5 = @(x) 1./(A.*( (3.975.*M) ./ (10.^6.* (x).^(3/2)) ).^n);

t5= integral(fun5, a05, af5);
