% ASE 324L Lab 6
% Anna Ring - aer3965

data = readmatrix('ASE324L_Lab6_Data.xlsx');

widths = data(1:6, 26);
thickness = data(1:6, 27);
widths = widths./1000; % m
thickness = thickness./1000; % m 
Area = widths.*thickness;
L = 50.8; % mm 
data = data(:, 1:24);

% Conversion factors
in_to_mm = 25.4;
lbf_to_N = 4.44822;

sample1 = data(:, 1:5);
sample2 = data(:, 6:10);
sample3 = data(:, 11:15);
sample4 = data(:, 16:18);
sample5 = data(:, 19:21);
sample6 = data(:, 22:24);

clear data

% Conversions
sample1(:,2) = sample1(:,2) * in_to_mm; % Extension: inches to mm
sample1(:,3) = sample1(:,3) * lbf_to_N; % Load: lbf to N

sample2(:,2) = sample2(:,2) * in_to_mm;
sample2(:,3) = sample2(:,3) * lbf_to_N;

sample3(:,2) = sample3(:,2) * in_to_mm;
sample3(:,3) = sample3(:,3) * lbf_to_N;

sample4(:,2) = sample4(:,2) * in_to_mm;
sample4(:,3) = sample4(:,3) * lbf_to_N;

sample5(:,2) = sample5(:,2) * in_to_mm;
sample5(:,3) = sample5(:,3) * lbf_to_N;

sample6(:,2) = sample6(:,2) * in_to_mm;
sample6(:,3) = sample6(:,3) * lbf_to_N;

% Stress, strain, etc
sample1(:, 3) = sample1(:, 3)./Area(1);
sample2(:, 3) = sample2(:, 3)./Area(2);
sample3(:, 3) = sample3(:, 3)./Area(3);
sample4(:, 3) = sample4(:, 3)./Area(4);
sample5(:, 3) = sample5(:, 3)./Area(5);
sample6(:, 3) = sample6(:, 3)./Area(6);

sample1(:, 2) = sample1(:, 2)./L;
sample2(:, 2) = sample2(:, 2)./L;
sample3(:, 2) = sample3(:, 2)./L;
sample4(:, 2) = sample4(:, 2)./L;
sample5(:, 2) = sample5(:, 2)./L;
sample6(:, 2) = sample6(:, 2)./L;


sample1 = truncateAtNaN(sample1);
sample2 = truncateAtNaN(sample2);
sample3 = truncateAtNaN(sample3);
sample4 = truncateAtNaN(sample4);
sample5 = truncateAtNaN(sample5);
sample6 = truncateAtNaN(sample6);

sample1 = remove_beginning_data(sample1);
sample2 = remove_beginning_data(sample2);
sample3 = remove_beginning_data(sample3);
sample4 = remove_beginning_data(sample4);
sample5 = remove_beginning_data(sample5);
sample6 = remove_beginning_data(sample6);

sample1(:, 2) = sample1(:, 2) - sample1(1,2);
sample2(:, 2) = sample2(:, 2) - sample2(1, 2);
sample3(:, 2) = sample3(:, 2) - sample3(1, 2);
sample4(:, 2) = sample4(:, 2) - sample4(1, 2);
sample5(:, 2) = sample5(:, 2) - sample5(1, 2);
sample6(:, 2) = sample6(:, 2) - sample6(1, 2);

sample1(:, 5) = sample1(:, 5) - sample1(1, 5);
sample2(:, 5) = sample2(:, 5) - sample2(1, 5);
sample3(:, 5) = sample3(:, 5) - sample3(1, 5);
sample1(:, 4) = sample1(:, 4) - sample1(1, 4);
sample2(:, 4) = sample2(:, 4) - sample2(1, 4);
sample3(:, 4) = sample3(:, 4) - sample3(1, 4);

% Find ultimate tensile strengths (maximum stress before failure)
UTS_0deg = max(sample4(:,3));  % 0 degree sample (parallel to fibers)
UTS_45deg = max(sample5(:,3)); % 45 degree sample
UTS_90deg = max(sample6(:,3)); % 90 degree sample (perpendicular to fibers)

% Display strengths
disp('--- Ultimate Strengths ---');
fprintf('Ultimate Tensile Strength (0° - parallel to fibers): %.2e Pa\n', UTS_0deg);
fprintf('Ultimate Tensile Strength (45°): %.2e Pa\n', UTS_45deg);
fprintf('Ultimate Tensile Strength (90° - perpendicular to fibers): %.2e Pa\n', UTS_90deg);


Enet1 = polyfit(sample1(:, 5), sample1(:, 3), 1);
Enet2 = polyfit(sample3(:, 5), sample3(:, 3), 1);
v21 = polyfit(sample1(:, 5), sample1(:, 4), 1);
v12 = polyfit(sample3(:, 5), sample3(:, 4), 1);
Exy = polyfit(sample2(:, 5), sample2(:, 3), 1);
Vxy = polyfit(sample2(:, 5), sample2(:, 4), 1);

v21 = -1*v21;
v12 = -1*v12;
Vxy = -1*Vxy;

% After fitting Young's Modulus and Poisson's Ratio

% Display results
disp('--- 0° Specimen (Sample 1) ---');
fprintf('Young''s Modulus (E1): %.2e Pa\n', Enet1(1));
fprintf('Poisson''s Ratio (v21): %.4f\n\n', v21(1));

disp('--- 45° Specimen (Sample 2) ---');
fprintf('Young''s Modulus (Exy): %.2e Pa\n', Exy(1));
fprintf('Poisson''s Ratio (vxy): %.4f\n\n', Vxy(1));
fprintf('Lamina Shear Modulus (G12): %.2e Pa\n\n', G(1));

disp('--- 90° Specimen (Sample 3) ---');
fprintf('Young''s Modulus (E2): %.2e Pa\n', Enet2(1));
fprintf('Poisson''s Ratio (v12): %.4f\n\n', v12(1));

disp('--- Volume Fraction ---');
fprintf('Fiber Volume Fraction (Vf): %.4f\n', volumeFraction);

G = Exy/(2*(1+Vxy(1)));

sMatrix = [1./Enet1(1), -v21(1)./Enet2(1), 0; -v12(1)./Enet1(1), 1./Enet2(1), 0; 0, 0, 1./G(1)];

volumeEqua = @(f)equaLab6(f, Enet1(1));
options = optimset('Display','off');
volumeFraction = fsolve(volumeEqua, 0.5, options);

%{
Ef = 235*10^9;
Em = 3*10^9;
Vm = 0.3;
Vf = 0.2;
f = 0.63;

v12Theory = (Vf*f*Em + (1-f)*Ef*Vm)/(f*Em + (1-f)*Ef);
%v12Theory = (Vf*f + Vm*(1-f))*(Ef*f + Em*(1-f))/inv((f/Ef) + ((1-f)/Em));
%}

figure;
hold on;
plot(sample1(:, 5), sample1(:, 3), 'LineWidth', 1.5);
plot(sample2(:, 5), sample2(:, 3), 'LineWidth', 1.5);
plot(sample3(:, 5), sample3(:, 3), 'LineWidth', 1.5);
xlabel('Strain')
ylabel('Engineering Stress (Pa)')
title('Stress Strain Relationship of Linear Regime')
legend('0 Deg Angle', '45 Deg Angle', '90 Deg Angle','Location','best');
hold off

figure;
hold on;
plot(sample4(:, 2), sample4(:, 3), 'LineWidth', 1.5);
plot(sample5(:, 2), sample5(:, 3), 'LineWidth', 1.5);
plot(sample6(:, 2), sample6(:, 3), 'LineWidth', 1.5);
xlabel('Strain')
ylabel('Engineering Stress (Pa)')
title('Stress Strain Relationship until Failure')
legend('0 Deg Angle', '45 Deg Angle', '90 Deg Angle','Location','best');
hold off


xVals1 = sample1(:, 5);
yVals1 = -1*(v21(1).*sample1(:, 5) + v21(2));

figure;
hold on;
plot(sample1(:, 5), sample1(:, 4), 'LineWidth', 1.5)
plot(xVals1, yVals1, 'LineWidth', 1.5)
xlabel('Axial Strain')
ylabel('Transverse Strain')
title('0 Degree Angle Sample Possions Ratio')
legend('Strain Data', 'Projected Possions Ratio', 'Location', 'best')
hold off;

xVals2 = sample2(:, 5);
yVals2 = -1*(Vxy(1).*sample2(:, 5) + Vxy(2));
figure;
hold on;
plot(sample2(:, 5), sample2(:, 4), 'LineWidth', 1.5)
plot(xVals2, yVals2, 'LineWidth', 1.5)
xlabel('Axial Strain')
ylabel('Transverse Strain')
title('45 Degree Angle Sample Possions Ratio')
legend('Strain Data', 'Projected Possions Ratio', 'Location', 'best')
hold off;

xVals3 = sample3(:, 5);
yVals3 = -1*(v12(1).*sample3(:, 5) + v12(2));
figure;
hold on;
plot(sample3(:, 5), sample3(:, 4), 'LineWidth', 1.5)
plot(xVals3, yVals3, 'LineWidth', 1.5)
xlabel('Axial Strain')
ylabel('Transverse Strain')
title('90 Degree Angle Sample Possions Ratio')
legend('Strain Data', 'Projected Possions Ratio', 'Location', 'best')
hold off;

figure;
hold on;
plot(sample1(:, 5), sample1(:, 3), 'LineWidth', 1.5)
plot(sample1(:, 5), Enet1(1).*sample1(:, 5) + Enet1(2), 'LineWidth', 1.5)
xlabel('Axial Strain')
ylabel('Engineering Stress (Pa)')
title('0 Degree Angle Sample Stress-Strain')
legend('Stress Data', 'Projected Youngs Modulus', 'Location', 'best')
hold off;

figure;
hold on;
plot(sample2(:, 5), sample2(:, 3), 'LineWidth', 1.5)
plot(sample2(:, 5), Exy(1).*sample2(:, 5) + Exy(2), 'LineWidth', 1.5)
xlabel('Axial Strain')
ylabel('Engineering Stress (Pa)')
title('45 Degree Angle Sample Stress-Strain')
legend('Stress Data', 'Projected Youngs Modulus', 'Location', 'best')
hold off;

figure;
hold on;
plot(sample3(:, 5), sample3(:, 3), 'LineWidth', 1.5)
plot(sample3(:, 5), Enet2(1).*sample3(:, 5) + Enet2(2), 'LineWidth', 1.5)
xlabel('Axial Strain')
ylabel('Engineering Stress (Pa)')
title('90 Degree Angle Sample Stress-Strain')
legend('Stress Data', 'Projected Youngs Modulus', 'Location', 'best')
hold off;

clear xVals1 xVals2 xVals3 yVals1 yVals2 yVals3

function truncatedMat = truncateAtNaN(mat)
    % Find the index of the first NaN in the vector
    vec = mat(:, 3);
    nanIndex = find(isnan(vec), 1);
    
    if isempty(nanIndex)
        % If no NaN is found, return the original vector
        truncatedMat = mat;
    else
        % Truncate the vector up to the first NaN
        truncatedMat = mat(1:nanIndex-1, :);
    end
end

function trimmed_array = remove_beginning_data(input_array)
    % Initialize index variable
    idx = 1;
    input_array_check = input_array(:, 3);
    % Iterate through the input array
    while idx <= length(input_array_check)
        % Check if current value is less than or equal to 0.01
        if input_array_check(idx) >= 0.03
            % If so, break out of the loop
            break;
        end
        % Increment index
        idx = idx + 1;
    end
    
    % Trim the input array from the beginning to the index where value <= 0.01
    trimmed_array = input_array(idx:end, :);
end