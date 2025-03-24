% Materials Lab 4 - Anna Ring
clc; close all; 

L = 33/1000; % m 
W = 12/1000; % m
T = 0.5/1000; % m
c = T/2;

excelData = readmatrix("Ceramics Data.xlsx");

Sample1Raw = excelData(1:(4958-10), 5:7);
Sample2Raw = excelData(1:(4594-10), 9:11);
Sample3Raw = excelData(1:(8574-10), 13:15);
Sample4Raw = excelData(1:(3772-10), 17:19);
Sample5Raw = excelData(1:(2930-10), 21:23);
Sample6Raw = excelData(1:(1299-10), 25:27);
Sample7Raw = excelData(1:(4719-10), 29:31);
Sample8Raw = excelData(1:(4191-10), 33:35);
Sample9Raw = excelData(1:(3072-10), 37:39);
Sample10Raw = excelData(1:(2071-10), 41:43);
Sample11Raw = excelData(1:(3587-10), 45:47);
Sample12Raw = excelData(1:(3768-10), 49:51);
Sample13Raw = excelData(1:(4485-10), 53:55);
Sample14Raw = excelData(1:(4122-10), 57:59);
Sample1 = remove_beginning_data(Sample1Raw);
Sample2 = remove_beginning_data(Sample2Raw);
Sample3 = remove_beginning_data(Sample3Raw);
Sample4 = remove_beginning_data(Sample4Raw);
Sample5 = remove_beginning_data(Sample5Raw);
Sample6 = remove_beginning_data(Sample6Raw);
Sample7 = remove_beginning_data(Sample7Raw);
Sample8 = remove_beginning_data(Sample8Raw);
Sample9 = remove_beginning_data(Sample9Raw);
Sample10 = remove_beginning_data(Sample10Raw);
Sample11 = remove_beginning_data(Sample11Raw);
Sample12 = remove_beginning_data(Sample12Raw);
Sample13 = remove_beginning_data(Sample13Raw);
Sample14 = remove_beginning_data(Sample14Raw);
Sample1(:, 2) = (12*c/L^2)*0.0254*(Sample1(:, 2) - Sample1(1, 2));
Sample2(:, 2) = (12*c/L^2)*0.0254*(Sample2(:, 2) - Sample2(1, 2));
Sample3(:, 2) = (12*c/L^2)*0.0254*(Sample3(:, 2) - Sample3(1, 2));
Sample4(:, 2) = (12*c/L^2)*0.0254*(Sample4(:, 2) - Sample4(1, 2));
Sample5(:, 2) = (12*c/L^2)*0.0254*(Sample5(:, 2) - Sample5(1, 2));
Sample6(:, 2) = (12*c/L^2)*0.0254*(Sample6(:, 2) - Sample6(1, 2));
Sample7(:, 2) = (12*c/L^2)*0.0254*(Sample7(:, 2) - Sample7(1, 2));
Sample8(:, 2) = (12*c/L^2)*0.0254*(Sample8(:, 2) - Sample8(1, 2));
Sample9(:, 2) = (12*c/L^2)*0.0254*(Sample9(:, 2) - Sample9(1, 2));
Sample10(:, 2) = (12*c/L^2)*0.0254*(Sample10(:, 2) - Sample10(1, 2));
Sample11(:, 2) = (12*c/L^2)*0.0254*(Sample11(:, 2) - Sample11(1, 2));
Sample12(:, 2) = (12*c/L^2)*0.0254*(Sample12(:, 2) - Sample12(1, 2));
Sample13(:, 2) = (12*c/L^2)*0.0254*(Sample13(:, 2) - Sample13(1, 2));
Sample14(:, 2) = (12*c/L^2)*0.0254*(Sample14(:, 2) - Sample14(1, 2));
Sample1(:, 3) = 4.44822*((3/2)*L/(W*T^2))*Sample1(:, 3);
Sample2(:, 3) = 4.44822*((3/2)*L/(W*T^2))*Sample2(:, 3);
Sample3(:, 3) = 4.44822*((3/2)*L/(W*T^2))*Sample3(:, 3);
Sample4(:, 3) = 4.44822*((3/2)*L/(W*T^2))*Sample4(:, 3);
Sample5(:, 3) = 4.44822*((3/2)*L/(W*T^2))*Sample5(:, 3);
Sample6(:, 3) = 4.44822*((3/2)*L/(W*T^2))*Sample6(:, 3);
Sample7(:, 3) = 4.44822*((3/2)*L/(W*T^2))*Sample7(:, 3);
Sample8(:, 3) = 4.44822*((3/2)*L/(W*T^2))*Sample8(:, 3);
Sample9(:, 3) = 4.44822*((3/2)*L/(W*T^2))*Sample9(:, 3);
Sample10(:, 3) = 4.44822*((3/2)*L/(W*T^2))*Sample10(:, 3);
Sample11(:, 3) = 4.44822*((3/2)*L/(W*T^2))*Sample11(:, 3);
Sample12(:, 3) = 4.44822*((3/2)*L/(W*T^2))*Sample12(:, 3);
Sample13(:, 3) = 4.44822*((3/2)*L/(W*T^2))*Sample13(:, 3);
Sample14(:, 3) = 4.44822*((3/2)*L/(W*T^2))*Sample14(:, 3);
Sample1E = polyfit(Sample1(:, 2), Sample1(:, 3), 1);
Sample2E = polyfit(Sample2(:, 2), Sample2(:, 3), 1);
Sample3E = polyfit(Sample3(:, 2), Sample3(:, 3), 1);
Sample4E = polyfit(Sample4(:, 2), Sample4(:, 3), 1);
Sample5E = polyfit(Sample5(:, 2), Sample5(:, 3), 1);
Sample6E = polyfit(Sample6(:, 2), Sample6(:, 3), 1);
Sample7E = polyfit(Sample7(:, 2), Sample7(:, 3), 1);
Sample8E = polyfit(Sample8(:, 2), Sample8(:, 3), 1);
Sample9E = polyfit(Sample9(:, 2), Sample9(:, 3), 1);
Sample10E = polyfit(Sample10(:, 2), Sample10(:, 3), 1);
Sample11E = polyfit(Sample11(:, 2), Sample11(:, 3), 1);
Sample12E = polyfit(Sample12(:, 2), Sample12(:, 3), 1);
Sample13E = polyfit(Sample13(:, 2), Sample13(:, 3), 1);
Sample14E = polyfit(Sample14(:, 2), Sample14(:, 3), 1);
SampleE = [Sample1E(1, 1), Sample2E(1, 1), Sample3E(1, 1), Sample4E(1, 1), Sample5E(1, 1), Sample6E(1, 1), Sample7E(1, 1), Sample8E(1, 1), Sample9E(1, 1), Sample10E(1, 1), Sample11E(1, 1), Sample12E(1, 1), Sample13E(1, 1), Sample14E(1, 1)];
SampleAvgE = mean(SampleE);

disp('Average E (Pa):')
disp(SampleAvgE)

Sample1FS = max(Sample1(:, 3));
Sample2FS = max(Sample2(:, 3));
Sample3FS = max(Sample3(:, 3));
Sample4FS = max(Sample4(:, 3));
Sample5FS = max(Sample5(:, 3));
Sample6FS = max(Sample6(:, 3));
Sample7FS = max(Sample7(:, 3));
Sample8FS = max(Sample8(:, 3));
Sample9FS = max(Sample9(:, 3));
Sample10FS = max(Sample10(:, 3));
Sample11FS = max(Sample11(:, 3));
Sample12FS = max(Sample12(:, 3));
Sample13FS = max(Sample13(:, 3));
Sample14FS = max(Sample14(:, 3));
SampleFS = [Sample1FS, Sample2FS, Sample3FS, Sample4FS, Sample5FS, Sample6FS, Sample7FS, Sample8FS, Sample9FS, Sample10FS, Sample11FS, Sample12FS, Sample13FS, Sample14FS];
SampleAvgFS = mean(SampleFS);
SampleFS = sort(SampleFS);

disp('Average FS (Pa):')
disp(SampleAvgFS)

PMeasured = 1:-(1/14):(1/14);

%stressValues = 0:20:max(SampleFS);

figure;
hold on;
plot(Sample1(:, 2), Sample1(:, 3))
plot(Sample2(:, 2), Sample2(:, 3))
plot(Sample3(:, 2), Sample3(:, 3))
plot(Sample4(:, 2), Sample4(:, 3))
plot(Sample5(:, 2), Sample5(:, 3))
plot(Sample6(:, 2), Sample6(:, 3))
plot(Sample7(:, 2), Sample7(:, 3))
plot(Sample8(:, 2), Sample8(:, 3))
plot(Sample9(:, 2), Sample9(:, 3))
plot(Sample10(:, 2), Sample10(:, 3))
plot(Sample11(:, 2), Sample11(:, 3))
plot(Sample12(:, 2), Sample12(:, 3))
plot(Sample13(:, 2), Sample13(:, 3))
plot(Sample14(:, 2), Sample14(:, 3))
xlabel('Strain from initial displacement')
ylabel('Stress (Pa)')
title('Stress Strain Diagrams')
legend('Sample 1', 'Sample 2','Sample 3','Sample 4','Sample 5','Sample 6','Sample 7','Sample 8' ...
    ,'Sample 9','Sample 10','Sample 11','Sample 12','Sample 13','Sample 14', 'Location', 'best')

% Define parameters
sigma_0 = 100; % Characteristic strength in MPa
m_values = [5, 10, 20]; % Weibull moduli

% Range of strength values
strength_values = linspace(0, 200, 1000); % Adjust the range as needed

% Plotting
figure;
hold on;

for i = 1:length(m_values)
    m = m_values(i);
    
    % Calculate Weibull distribution
    %weibull_pdf = (m / sigma_0) * (strength_values / sigma_0).^(m - 1) .* exp(-(strength_values / sigma_0).^m);
    weibull_pdf = exp(-(strength_values./sigma_0).^m);
    % Plot Weibull distribution
    plot(strength_values, weibull_pdf, 'LineWidth', 2, 'DisplayName', ['m = ', num2str(m)], 'Color', rand(3, 1));
end

hold off;

% Add labels and legend
xlabel('Applied Strength (Pa)');
ylabel('Probability of survival');
title('Weibull Distribution of Ceramic Strength');
legend('Location', 'best');

failure_data = sort(SampleFS);
n = length(failure_data);

survival_prob = (n - (0.5:(n-0.5)))/n;
y_vals = log(log(1./survival_prob));
x_vals = log(failure_data);
coeff = polyfit(x_vals, y_vals, 1);
m = coeff(1);
sigma_0_fit = exp(-coeff(2)/m);

stressVals = linspace(0, 2*sigma_0_fit, 1000);
%stressVals = linspace(1.5e+8, 4.5e+8, 1000);
survival_Prob_fit = exp(-(stressVals./sigma_0_fit).^m);

figure;
hold on;
plot(stressVals, survival_Prob_fit, 'LineWidth', 1.5)
plot(failure_data, survival_prob, '*')
xlabel('Stress (Pa)')
ylabel('Survival Probablity')
title('Weibull Distribution of Measured Samples')
legend('Projected Survival', 'Actual Data', 'Location', 'best')
hold off;

% Checking question 5
R5 = 5/1000;
L5 = 15/1000;
d5 = 0.01/1000;
F5 = 7500;
strain5 = 12*R5*d5/(L5*L5);
stress5 = F5*L5/(pi*R5^3);
E5 = stress5/strain5;


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