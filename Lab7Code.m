% ASE 324L Lab 7 Code
% Anna Ring - aer3965

data = readmatrix("Lab_7_fract.xlsx");
inToM = 0.0254;
lbfToN = 4.44822;
W = 2;
W = W*inToM;

b = [0.5, 0.25, 0.125];
b = b.*inToM;

a05 = data(1, 2:2:12);
a025 = data(1, 15:2:23);
a0125 = data(1, 26:2:36);
a = [a05', [a025, 0]', a0125'];
a = a.*inToM;

data = data(4:end, :);
data05 = data(:, 1:12);
data025 = data(:, 14:23);
data0125 = data(:, 25:end);

clear data a05 a025 a0125

for ii = 1:6
    data05(:, 2*ii) = lbfToN.*data05(:, 2*ii);
    data05(:, 2*ii-1) = inToM.*data05(:, 2*ii-1);
    data0125(:, 2*ii) = lbfToN.*data0125(:, 2*ii);
    data0125(:, 2*ii-1) = inToM.*data0125(:, 2*ii-1);

    if ii > 5.9
        break;
    else
        data025(:, 2*ii) = lbfToN.*data025(:, 2*ii);
        data025(:, 2*ii-1) = inToM.*data025(:, 2*ii-1);
    end
end
Kc = 0.*a;
x = a./W;
aVec = 29.6.*sqrt(x) - 185.5.*(x).^(3/2) + 655.7.*(x).^(5/2) - 1017.*(x).^(7/2) + 639.*(x).^(9/2);

for ii = 1:length(a)
    Kc(ii, 1) = max(data05(:, 2*ii)).*(aVec(ii, 1))./(b(1).*sqrt(W));
    Kc(ii, 3) = max(data0125(:, 2*ii)).*(aVec(ii, 3))./(b(3).*sqrt(W));
    if ii > 5.9
        Kc(ii, 2) = NaN;
        break;
    else
        Kc(ii, 2) = max(data025(:, 2*ii)).*(aVec(ii, 2))./(b(2).*sqrt(W));
    end
end
a(end, 2) = NaN;

% Convert Kc from Pa to MPa for display
Kc_MPa = Kc ./ 1e6;

% Create row names and column labels
sampleLabels = ["Sample 1"; "Sample 2"; "Sample 3"; "Sample 4"; "Sample 5"; "Sample 6"];
colLabels = ["B = 0.5 in", "B = 0.25 in", "B = 0.125 in"];

% Create tables for a values and Kc values
a_table = array2table(a, 'VariableNames', colLabels, 'RowNames', sampleLabels);
Kc_table = array2table(Kc_MPa, 'VariableNames', colLabels, 'RowNames', sampleLabels);

% Display the tables
disp('A Values (m):');
disp(a_table);

disp('Kc Values (MPa):');
disp(Kc_table);


ysAl = 57; %ksi
ksiToPa = 6.895e+6;
ysAl = ysAl*ksiToPa;

KicMetric = 2.5.*(Kc./ysAl).^2;
KicMetricBool = KicMetric <= a;

KicVals = KicMetricBool.*Kc;


inputHalf = zeros(length(data0125), 2);

C05 = zeros(6, 1);
C025 = zeros(5, 1);
C125 = zeros(6, 1);
CInter = 0;

for ii = 1:6
    inputHalf = data05(:, (2*ii-1):(2*ii));
    inputHalf = remove_beginning_data(inputHalf);
    inputHalf(:, 1) = inputHalf(:, 1) - inputHalf(1, 1);
    data05(:, (2*ii-1):(2*ii)) = [inputHalf; NaN(length(data05) - length(inputHalf) , 2)];
    inputHalf = truncateAtNaN(inputHalf);
    CInter = polyfit(inputHalf(200:500, 2), inputHalf(200:500, 1), 1);
    C05(ii) = CInter(1);

    
    if ii < 5.9
        inputHalf = data025(:, (2*ii-1):(2*ii));
        inputHalf = remove_beginning_data(inputHalf);
        inputHalf(:, 1) = inputHalf(:, 1) - inputHalf(1, 1);
        data025(:, (2*ii-1):(2*ii)) = [inputHalf; NaN(length(data05) - length(inputHalf) , 2)];
        inputHalf = truncateAtNaN(inputHalf);
        CInter = polyfit(inputHalf(200:500, 2), inputHalf(200:500, 1), 1);
        C025(ii) = CInter(1);
    end

    inputHalf = data0125(:, (2*ii-1):(2*ii));
    inputHalf = remove_beginning_data(inputHalf);
    inputHalf(:, 1) = inputHalf(:, 1) - inputHalf(1, 1);
    data0125(:, (2*ii-1):(2*ii)) = [inputHalf; NaN(length(data05) - length(inputHalf) , 2)];
    inputHalf = truncateAtNaN(inputHalf);
    CInter = polyfit(inputHalf(200:500, 2), inputHalf(200:500, 1), 1);
    C125(ii) = CInter(1);
end

C = [C05, [C025; NaN], C125];

% Sort the "a" matrix
a_sorted = a.*0;
Kc_sorted = Kc.*0;
Kc_Mean = [0, 0, 0];
C_sorted = C.*0;
for ii = 1:3
    [a_sorted(:, ii), idx] = sortrows(a(:, ii));
    % Sort the "Kc" matrix by the same changes
    Kc_sorted(:, ii) = Kc(idx, ii);
    C_sorted(:, ii) = C(idx, ii);
    if ii < 2.1 && ii > 1.9
        Kc_Mean(ii) = mean(Kc(1:5, ii));
    else
        Kc_Mean(ii) = mean(Kc(:, ii));
    end
end


disp('Average Kc (Pa) for each thickness:');
disp(array2table(Kc_Mean, 'VariableNames', {'B_0_5_in', 'B_0_25_in', 'B_0_125_in'}));


figure;
hold on;
plot(a_sorted(:, 1), Kc_sorted(:, 1), 'LineWidth', 1.5)
plot(a_sorted(1:5, 2), Kc_sorted(1:5, 2), 'LineWidth', 1.5)
plot(a_sorted(:, 3), Kc_sorted(:, 3), 'LineWidth', 1.5)
xlabel('Crack Length (m)')
ylabel('Kc')
title('Kc v a')
legend('0.5 in Thickness', '0.25 in Thickness', '0.125 in Thickness', 'Location','best')
hold off;

figure;
hold on;
plot(b, Kc_Mean, 'LineWidth', 1.5)
xlabel("Thickness (m)")
ylabel("Average Kc (Pa)")
title('Average Toughness')
hold off;

figure;
hold on;
for ii = 1:6
    hold on;
    plot(data05(:, 2*ii-1), data05(:, 2*ii), 'LineWidth', 1.5)
end
xlabel('Displacement (m)')
ylabel('Forcing Load (N)')
title('Load by Displacement for 0.5 in Sample ')
legend('Sample 1','Sample 2','Sample 3','Sample 4','Sample 5','Sample 6', 'Location', 'best')
hold off

figure;
hold on;
for ii = 1:5
    hold on;
    plot(data025(:, 2*ii-1), data025(:, 2*ii), 'LineWidth', 1.5)
end
xlabel('Displacement (m)')
ylabel('Forcing Load (N)')
title('Load by Displacement for 0.25 in Sample ')
legend('Sample 1','Sample 2','Sample 3','Sample 4','Sample 5', 'Location', 'best')
hold off

figure;
hold on;
for ii = 1:6
    hold on;
    plot(data0125(:, 2*ii-1), data0125(:, 2*ii), 'LineWidth', 1.5)
end
xlabel('Displacement (m)')
ylabel('Forcing Load (N)')
title('Load by Displacement for 0.125 in Sample ')
legend('Sample 1','Sample 2','Sample 3','Sample 4','Sample 5','Sample 6', 'Location', 'best')
hold off

figure;
hold on;
plot(a_sorted(:, 1), C_sorted(:, 1), 'LineWidth', 1.5)
plot(a_sorted(:, 2), C_sorted(:, 2), 'LineWidth', 1.5)
plot(a_sorted(:, 3), C_sorted(:, 3), 'LineWidth', 1.5)
xlabel('Crack Length (m)')
ylabel('Compliance (Pa)^-1')
title('Compliance Over Crack Thickness')
legend('0.5 in Thickness', '0.25 in Thickness', '0.125 in Thickness', 'Location','best')
hold off;

function trimmed_array = remove_beginning_data(input_array)
    % Initialize index variable
    idx = 1;
    input_array_check = input_array(:, 2);
    % Iterate through the input array
    while idx <= length(input_array_check)
        % Check if current value is less than or equal to 0.01
        if input_array_check(idx) >= 50
            % If so, break out of the loop
            break;
        end
        % Increment index
        idx = idx + 1;
    end
    
    % Trim the input array from the beginning to the index where value <= 0.01
    trimmed_array = input_array(idx:end, :);
end

function truncatedMat = truncateAtNaN(mat)
    % Find the index of the first NaN in the vector
    vec = mat(:, 2);
    nanIndex = find(isnan(vec), 1);
    
    if isempty(nanIndex)
        % If no NaN is found, return the original vector
        truncatedMat = mat;
    else
        % Truncate the vector up to the first NaN
        truncatedMat = mat(1:nanIndex-1, :);
    end
end