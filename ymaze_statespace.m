clc;
close all;

% Load the first dataset
filename1 = 'C:\Users\Anacker1\Desktop\logistic regression\vCA1-R-INH-Cont.xlsx'; % Update with your actual file path
opts1 = detectImportOptions(filename1);
opts1.DataRange = 'A2';  % Start importing data from the second row
data1 = readmatrix(filename1, opts1);
trials1 = data1(:, 1);
responses1 = data1(:, 2:end);

% Initialize matrix for storing predicted probabilities using state-space model
numTrials1 = length(trials1);
numMice1 = size(responses1, 2);
predictedProbabilitiesSSM1 = zeros(numTrials1, numMice1);

% State-space model prediction for each mouse in the first dataset
for i = 1:numMice1
    A = 1; % State transition matrix
    C = 1; % Observation matrix
    Q = 0.01; % Process noise covariance
    R = 0.1; % Measurement noise covariance                                                            
    x = 0.5; % Initial state estimate
    P = 1; % Initial error covariance

    % Kalman filter implementation
    for t = 1:numTrials1
        x = A * x; % Predicted state estimate
        P = A * P * A' + Q; % Predicted error covariance
        K = P * C' / (C * P * C' + R); % Kalman gain
        x = x + K * (responses1(t, i) - C * x); % Updated state estimate
        P = (1 - K * C) * P; % Updated error covariance
        predictedProbabilitiesSSM1(t, i) = x;
    end
end

% Aggregate the predictions across mice for each trial
meanProbabilitiesSSM1 = mean(predictedProbabilitiesSSM1, 2);
windowSize = 5; % Choose an appropriate window size based on your data
smoothedProbabilitiesSSM1 = movmean(meanProbabilitiesSSM1, windowSize);
SEM_SSM1 = std(predictedProbabilitiesSSM1, 0, 2) / sqrt(numMice1);
smoothedSEM_SSM1 = movmean(SEM_SSM1, windowSize);

trialRange1 = trials1;
upperBound_SSM1 = smoothedProbabilitiesSSM1 + smoothedSEM_SSM1;
lowerBound_SSM1 = smoothedProbabilitiesSSM1 - smoothedSEM_SSM1;

% Load the second dataset
filename2 = 'C:\Users\Anacker1\Desktop\logistic regression\vCA1-R-INH-Dreadd.xlsx'; % Update with your actual file path
opts2 = detectImportOptions(filename2);
opts2.DataRange = 'A2';  % Start importing data from the second row
data2 = readmatrix(filename2, opts2);
trials2 = data2(:, 1);
responses2 = data2(:, 2:end);

numTrials2 = length(trials2);
numMice2 = size(responses2, 2);
predictedProbabilitiesSSM2 = zeros(numTrials2, numMice2);

% State-space model prediction for each mouse in the second dataset
for i = 1:numMice2
    A = 1;
    C = 1;
    Q = 0.01;
    R = 0.1;
    x = 0.5;
    P = 1;

    for t = 1:numTrials2
        x = A * x;
        P = A * P * A' + Q;
        K = P * C' / (C * P * C' + R);
        x = x + K * (responses2(t, i) - C * x);
        P = (1 - K * C) * P;
        predictedProbabilitiesSSM2(t, i) = x;
    end
end

meanProbabilitiesSSM2 = mean(predictedProbabilitiesSSM2, 2);
smoothedProbabilitiesSSM2 = movmean(meanProbabilitiesSSM2, windowSize);
SEM_SSM2 = std(predictedProbabilitiesSSM2, 0, 2) / sqrt(numMice2);
smoothedSEM_SSM2 = movmean(SEM_SSM2, windowSize);

trialRange2 = trials2;
upperBound_SSM2 = smoothedProbabilitiesSSM2 + smoothedSEM_SSM2;
lowerBound_SSM2 = smoothedProbabilitiesSSM2 - smoothedSEM_SSM2;

% Combine both ranges (1-35 and 35-end) in a single figure
figure;

% Plot trials 1 to 35 for both datasets
range1_1to35 = trialRange1 >= 1 & trialRange1 <= 35;
range2_1to35 = trialRange2 >= 1 & trialRange2 <= 35;

fill([trialRange1(range1_1to35); flipud(trialRange1(range1_1to35))], [upperBound_SSM1(range1_1to35); flipud(lowerBound_SSM1(range1_1to35))], [0.9 0.9 1], 'linestyle', 'none');
hold on;
plot(trialRange1(range1_1to35), smoothedProbabilitiesSSM1(range1_1to35), 'b-', 'LineWidth', 2);

fill([trialRange2(range2_1to35); flipud(trialRange2(range2_1to35))], [upperBound_SSM2(range2_1to35); flipud(lowerBound_SSM2(range2_1to35))], [1 0.9 0.9], 'linestyle', 'none');
plot(trialRange2(range2_1to35), smoothedProbabilitiesSSM2(range2_1to35), 'r-', 'LineWidth', 2);

% Plot trials 35 to the end for both datasets
range1_35toEnd = trialRange1 >= 35;
range2_35toEnd = trialRange2 >= 35;

fill([trialRange1(range1_35toEnd); flipud(trialRange1(range1_35toEnd))], [upperBound_SSM1(range1_35toEnd); flipud(lowerBound_SSM1(range1_35toEnd))], [0.9 0.9 1], 'linestyle', 'none');
plot(trialRange1(range1_35toEnd), smoothedProbabilitiesSSM1(range1_35toEnd), 'b--', 'LineWidth', 2);

fill([trialRange2(range2_35toEnd); flipud(trialRange2(range2_35toEnd))], [upperBound_SSM2(range2_35toEnd); flipud(lowerBound_SSM2(range2_35toEnd))], [1 0.9 0.9], 'linestyle', 'none');
plot(trialRange2(range2_35toEnd), smoothedProbabilitiesSSM2(range2_35toEnd), 'r--', 'LineWidth', 2);

% Add labels and legend
xlabel('Trial Number');
ylabel('Probability of Correct Response');
title('Learning Curves with SEM Shading (State-Space Model)');
legend('First SEM (1-35)', 'First Learning Curve (1-35)', 'Second SEM (1-35)', 'Second Learning Curve (1-35)', ...
       'First SEM (35-End)', 'First Learning Curve (35-End)', 'Second SEM (35-End)', 'Second Learning Curve (35-End)');
grid on;
hold off;
p_values = zeros(numTrials1, 1); % Assuming numTrials1 == numTrials2 for simplicity
for t = 1:numTrials1
    [~, p_values(t)] = ttest2(predictedProbabilitiesSSM1(t, :), predictedProbabilitiesSSM2(t, :));
end

