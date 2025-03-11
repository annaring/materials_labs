clear all; close all; clc;
% ASE 324L Lab 3 Sp24 - Anna Ring

% I put all data in a separate csv file
un = readmatrix("untreated.csv");
an = readmatrix("annealed.csv");
a05 = readmatrix("aged0.5h.csv");
a2 = readmatrix('aged2h.csv');
a6 = readmatrix('aged6h.csv');
a24 = readmatrix('aged24h.csv');

unArea = (pi/4)*(0.502)^2;
anArea = (pi/4)*(0.5025)^2;
a05Area = (pi/4)*(0.499)^2;
a2Area = (pi/4)*(0.501)^2;
a6Area = (pi/4)*(0.501)^2;
a24Area = (pi/4)*(0.5025)^2;

L = 2; % inch 

% Strain 
un(:, 2) = un(:, 2)/L;
an(:, 2) = an(:, 2)/L;
a05(:, 2) = a05(:, 2)/L;
a2(:, 2) = a2(:, 2)/L;
a6(:, 2) = a6(:, 2)/L;
a24(:, 2) = a24(:, 2)/L;

% Stress
un(:, 1) = un(:, 1)/unArea;
an(:, 1) = an(:, 1)/anArea;
a05(:, 1) = a05(:, 1)/a05Area;
a2(:, 1) = a2(:, 1)/a2Area;
a6(:, 1) = a6(:, 1)/a6Area;
a24(:, 1) = a24(:, 1)/a24Area;

% UTS
unUTS = max(un(:, 1));
anUTS = max(an(:, 1));
a05UTS = max(a05(:, 1));
a2UTS = max(a2(:, 1));
a6UTS = max(a6(:, 1));
a24UTS = max(a24(:, 1));

% Toughness
unT = trapz(un(:, 3), un(:, 1)); 
anT = trapz(an(:, 3), an(:, 1));
a05T = trapz(a05(:, 3), a05(:, 1));
a2T = trapz(a2(:, 3), a2(:, 1));
a6T = trapz(a6(:, 3), a6(:, 1));
a24T = trapz(a24(:, 3), a24(:, 1));

unDA = (pi*(0.502/2)^2 - pi*(0.4/2)^2)/(pi*(0.502/2)^2);
unDL = (2.5065-2)/2;

anDA = (pi*(0.5025/2)^2 - pi*(0.4105/2)^2)/(pi*(0.5025/2)^2);
anDL = (2.2985-2)/2;

a05DA = (pi*(0.499/2)^2 - pi*(0.4035/2)^2)/(pi*(0.499/2)^2);
a05DL = (2.5275-2)/2;

a2DA = (pi*(0.501/2)^2 - pi*(0.401/2)^2)/(pi*(0.501/2)^2);
a2DL = (2.43-2)/2;

a6DA = (pi*(0.501/2)^2 - pi*(0.401/2)^2)/(pi*(0.501/2)^2);
a6DL = (2.495-2)/2;

a24DA = (pi*(0.5025/2)^2 - pi*(0.423/2)^2)/(pi*(0.5025/2)^2);
a24DL = (2.263-2)/2;

unVals = 300:800;
unSa = un(unVals, 3);
unSe = un(unVals, 1);
unET = -un(unVals, 4);
unE = polyfit(unSa, unSe, 1); 
unV = polyfit(unSa, unET, 1);

anVals = 400:800;
anSa = an(anVals, 3);
anSe = an(anVals, 1);
anET = -an(anVals, 4);
anE = polyfit(anSa, anSe, 1);
anV = polyfit(anSa, anET, 1);

a05Vals = 300:800;
a05Sa = a05(a05Vals, 3);
a05Se = a05(a05Vals, 1);
a05ET = -a05(a05Vals, 4);
a05E = polyfit(a05Sa, a05Se, 1);
a05V = polyfit(a05Sa, a05ET, 1);

a2Vals = 400:800;
a2Sa = a2(a2Vals, 3);
a2Se = a2(a2Vals, 1);
a2ET = -a2(a2Vals, 4);
a2E = polyfit(a2Sa, a2Se, 1);
a2V = polyfit(a2Sa, a2ET, 1);

a6Vals = 400:800;
a6Sa = a6(a6Vals, 3);
a6Se = a6(a6Vals, 1);
a6ET = -a6(a6Vals, 4);
a6E = polyfit(a6Sa, a6Se, 1);
a6V = polyfit(a6Sa, a6ET, 1);

a24Vals = 600:1000;
a24Sa = a24(a24Vals, 3);
a24Se = a24(a24Vals, 1);
a24ET = -a24(a24Vals, 4);
a24E = polyfit(a24Sa, a24Se, 1 );
a24V = polyfit(a24Sa, a24ET, 1);

% Display Ductility (Reduction in Area)
disp('Ductility (Reduction in Area):')
fprintf('Untreated: %.4f\n', unDA)
fprintf('Annealed: %.4f\n', anDA)
fprintf('Aged 0.5 hr: %.4f\n', a05DA)
fprintf('Aged 2 hr: %.4f\n', a2DA)
fprintf('Aged 6 hr: %.4f\n', a6DA)
fprintf('Aged 24 hr: %.4f\n', a24DA)

% Display Ductility (Elongation)
disp('Ductility (Elongation):')
fprintf('Untreated: %.4f\n', unDL)
fprintf('Annealed: %.4f\n', anDL)
fprintf('Aged 0.5 hr: %.4f\n', a05DL)
fprintf('Aged 2 hr: %.4f\n', a2DL)
fprintf('Aged 6 hr: %.4f\n', a6DL)
fprintf('Aged 24 hr: %.4f\n', a24DL)


%%
figure;
hold on;
plot(a24(:, 2), a24(:, 1), 'LineWidth', 1.5)
plot(a24(:, 3), a24(:, 1), 'LineWidth', 1.5)
title('Aged 24 Hours Aluminum Stress-Strain')
xlabel('Strain')
ylabel('Stress (psi)')
legend('Crosshead Displacement Strain', 'Extensometer Strain', 'Location', 'best')
hold off;

figure;
hold on;
plot(a6(:, 2), a6(:, 1), 'LineWidth', 1.5)
plot(a6(:, 3), a6(:, 1), 'LineWidth', 1.5)
title('Aged 6 Hours Aluminum Stress-Strain')
xlabel('Strain')
ylabel('Stress (psi)')
legend('Crosshead Displacement Strain', 'Extensometer Strain', 'Location', 'best')
hold off;

figure;
hold on;
plot(a2(:, 2), a2(:, 1), 'LineWidth', 1.5)
plot(a2(:, 3), a2(:, 1), 'LineWidth', 1.5)
title('Aged 2 Hours Aluminum Stress-Strain')
xlabel('Strain')
ylabel('Stress (psi)')
legend('Crosshead Displacement Strain', 'Extensometer Strain', 'Location', 'best')
hold off;

figure;
hold on;
plot(a05(:, 2), a05(:, 1), 'LineWidth', 1.5)
plot(a05(:, 3), a05(:, 1), 'LineWidth', 1.5)
title('Aged 0.5 Hours Aluminum Stress-Strain')
xlabel('Strain')
ylabel('Stress (psi)')
legend('Crosshead Displacement Strain', 'Extensometer Strain', 'Location', 'best')
hold off;

figure;
hold on;
plot(an(:, 2), an(:, 1), 'LineWidth', 1.5)
plot(an(:, 3), an(:, 1), 'LineWidth', 1.5)
title('Annealed Aluminum Stress-Strain')
xlabel('Strain')
ylabel('Stress (psi)')
legend('Crosshead Displacement Strain', 'Extensometer Strain', 'Location', 'best')
hold off;

figure;
hold on;
plot(un(:, 2), un(:, 1), 'LineWidth', 1.5)
plot(un(:, 3), un(:, 1), 'LineWidth', 1.5)
title('Untreated Aluminum Stress-Strain')
xlabel('Strain')
ylabel('Stress (psi)')
legend('Crosshead Displacement Strain', 'Extensometer Strain', 'Location', 'best')
hold off;

x = 0:0.0001:0.007;
figure;
hold on;
plot(un(:, 3), un(:, 1), 'LineWidth', 1.5)
plot(x, unE(1).*(x -0.000) + unE(2), 'LineWidth', 1.5)
title('Untreated Aluminum Bulk Modulus')
xlabel('Strain')
ylabel('Stress (psi)')
legend('Extensonmeter Strain', "Young's Modulus", 'Location', 'best')
hold off;

figure;
hold on;
plot(an(:, 3), an(:, 1), 'LineWidth', 1.5)
plot(x, anE(1).*(x -0.000) + anE(2), 'LineWidth', 1.5)
title('Annealed Aluminum Bulk Modulus')
xlabel('Strain')
ylabel('Stress (psi)')
legend('Extensometer Strain', "Young's Modulus", 'Location', 'best')
hold off;

figure;
hold on;
plot(a05(:, 3), a05(:, 1), 'LineWidth', 1.5)
plot(x, a05E(1).*(x -0.000) + a05E(2), 'LineWidth', 1.5)
title('Aged 0.5 Hours Aluminum Bulk Modulus')
xlabel('Strain')
ylabel('Stress (psi)')
legend('Extensometer Strain', "Young's Modulus", 'Location', 'best')
hold off;

figure;
hold on;
plot(a2(:, 3), a2(:, 1), 'LineWidth', 1.5)
plot(x, a2E(1).*(x -0.000) + a2E(2), 'LineWidth', 1.5)
title('Aged 2 Hours Aluminum Bulk Modulus')
xlabel('Strain')
ylabel('Stress (psi)')
legend('Extensometer Strain', "Young's Modulus", 'Location', 'best')
hold off;

figure;
hold on;
plot(a6(:, 3), a6(:, 1), 'LineWidth', 1.5)
plot(x, a6E(1).*(x -0.000) + a6E(2), 'LineWidth', 1.5)
title('Aged 6 Hours Aluminum Bulk Modulus')
xlabel('Strain')
ylabel('Stress (psi)')
legend('Extensometer Strain', "Young's Modulus", 'Location', 'best')
hold off;

figure;
hold on;
plot(a24(:, 3), a24(:, 1), 'LineWidth', 1.5)
plot(x, a24E(1).*(x -0.000) + a24E(2), 'LineWidth', 1.5)
title('Aged 24 Hours Aluminum Bulk Modulus')
xlabel('Strain')
ylabel('Stress (psi)')
legend('Extensometer Strain', "Young's Modulus", 'Location', 'best')
hold off;

xPos = 0:0.0001:0.005;
PosVals = 1:1000;
figure;
hold on;
plot(un(PosVals, 3), -un(PosVals, 4), 'LineWidth', 1.5)
plot(xPos, unV(1).*(xPos -0.000) + unV(2), 'LineWidth', 1.5)
title('Untreated Aluminum Poisson Ratio')
xlabel('Axial Strain')
ylabel('Transverse Strain')
legend('Strain Ratio', "Poisson Ratio", 'Location', 'best')
hold off;

figure;
hold on;
plot(an(PosVals, 3), -an(PosVals, 4), 'LineWidth', 1.5)
plot(xPos, anV(1).*(xPos -0.000) + anV(2), 'LineWidth', 1.5)
title('Annealed Aluminum Poisson Ratio')
xlabel('Axial Strain')
ylabel('Transverse Strain')
legend('Strain Ratio', "Poisson Ratio", 'Location', 'best')
hold off;

figure;
hold on;
plot(a05(PosVals, 3), -a05(PosVals, 4), 'LineWidth', 1.5)
plot(xPos, a05V(1).*(xPos -0.000) + a05V(2), 'LineWidth', 1.5)
title('Aged 0.5 Hours Aluminum Poisson Ratio')
xlabel('Axial Strain')
ylabel('Transverse Strain')
legend('Strain Ratio', "Poisson Ratio", 'Location', 'best')
hold off;

figure;
hold on;
plot(a2(PosVals, 3), -a2(PosVals, 4), 'LineWidth', 1.5)
plot(xPos, a2V(1).*(xPos -0.000) + a2V(2), 'LineWidth', 1.5)
title('Aged 2 Hours Aluminum Poisson Ratio')
xlabel('Axial Strain')
ylabel('Transverse Strain')
legend('Strain Ratio', "Poisson Ratio", 'Location', 'best')
hold off;

figure;
hold on;
plot(a6(PosVals, 3), -a6(PosVals, 4), 'LineWidth', 1.5)
plot(xPos, a6V(1).*(xPos -0.000) + a6V(2), 'LineWidth', 1.5)
title('Aged 6 Hours Aluminum Poisson Ratio')
xlabel('Axial Strain')
ylabel('Transverse Strain')
legend('Strain Ratio', "Poisson Ratio", 'Location', 'best')
hold off;

figure;
hold on;
plot(a24(PosVals, 3), -a24(PosVals, 4), 'LineWidth', 1.5)
plot(xPos, a24V(1).*(xPos -0.000) + a24V(2), 'LineWidth', 1.5)
title('Aged 24 Hours Aluminum Poisson Ratio')
xlabel('Axial Strain')
ylabel('Transverse Strain')
legend('Strain Ratio', "Poisson Ratio", 'Location', 'best')
hold off;

figure;
hold on;
plot(un(:, 3), un(:, 1), 'LineWidth', 1.5)
plot(an(:, 3), an(:, 1), 'LineWidth', 1.5)
plot(a05(:, 3), a05(:, 1), 'LineWidth', 1.5)
plot(a2(:, 3), a2(:, 1), 'LineWidth', 1.5)
plot(a6(:, 3), a6(:, 1), 'LineWidth', 1.5)
plot(a24(:, 3), a24(:, 1), 'LineWidth', 1.5)
legend('Untreated', 'Annealed', 'Aged 0.5 hr', 'Aged 2 hr', 'Aged 6 Hr', 'Aged 24 Hr', 'Location', 'best')
title('All Stress-Strain Curves')
xlabel('Strain (Extensometer)')
ylabel('Stress (psi)')
hold off;