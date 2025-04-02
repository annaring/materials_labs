% ASE324L Lab 5 - Polymers
% Anna Ring - aer3965

% Conversion from US to metric
inTom = 0.0254; 
lbfToN = 4.44822;

W = inTom*0.497;
T = inTom*0.256;
L = inTom*2;
A = W*T;

dataVals = readmatrix("High_density_PE.csv");
Temp = dataVals(7:16 , 1);
dataVals = dataVals(:, 2:31);
Er = zeros(length(Temp), 1);
Et = zeros(length(dataVals), 1);
time = dataVals(:, 1);

for ii = 1:10
    dataVals(:, 3*ii-1) = (inTom/L).*dataVals(:, 3*ii-1);
    dataVals(:, 3*ii) = (lbfToN/A).*dataVals(:, 3*ii);
    Er(ii) = dataVals(505, 3*ii)/dataVals(505, 3*ii -1);
end

disp('Relaxation modulus Er at t = 10s for each temperature:')
for i = 1:length(Temp)
    fprintf('Temperature = %.2f °C, Er = %.3e Pa\n', Temp(i), Er(i));
end

for ii = 1:length(dataVals)
    Et(ii) = dataVals(ii, 3)./dataVals(ii, 2);
end

Et = Et(2:6090, 1);

yVal = log(Et);
xVal = time(2:6090);

coeff = polyfit(xVal, yVal, 1);

E0 = exp(coeff(2) - 1);
eta = -E0/coeff(1);

fprintf('Spring constant E = %.3e Pa\n', E0);
fprintf('Viscosity η = %.3e Pa·s\n', eta);


figure;
plot(Temp, log(Er))
xlabel('Temperature (C)')
ylabel('log(Er) [log(Pa)]')
title('log of Relaxation Modulus at t = 10 s vs Temperature')