
clc;
close all;

% Load the first dataset
filename1 = 'C:\Users\Anacker1\Desktop\logistic regression\dCA1-R-INH-Cont.xlsx'; % Update with your actual file path
opts1 = detectImportOptions(filename1);
opts1.DataRange = 'A2';  % Start importing data from the second row
data1 = readmatrix(filename1, opts1);
trials1 = data1(:, 1);
responses1 = data1(:, 2:end);

% Load the second dataset
filename2 = 'C:\Users\Anacker1\Desktop\logistic regression\dCA1-R-INH-Dreadd.xlsx'; % Update with your actual file path
opts2 = detectImportOptions(filename2);
opts2.DataRange = 'A2';  % Start importing data from the second row
data2 = readmatrix(filename2, opts2);
trials2 = data2(:, 1);
responses2 = data2(:, 2:end);

% Analyze the first dataset for trials 1-36
range1_1_36 = trials1 >= 1 & trials1 <= 36;
[inflectionTrial1_1_36, meanProbabilities1_1_36, smoothedProbabilities1_1_36, upperBound1_1_36, lowerBound1_1_36] = analyzeLearningCurve(trials1(range1_1_36), responses1(range1_1_36, :));

% Analyze the first dataset for trials 36-end
range1_36_end = trials1 > 36;
[inflectionTrial1_36_end, meanProbabilities1_36_end, smoothedProbabilities1_36_end, upperBound1_36_end, lowerBound1_36_end] = analyzeLearningCurve(trials1(range1_36_end), responses1(range1_36_end, :));

% Analyze the second dataset for trials 1-36
range2_1_36 = trials2 >= 1 & trials2 <= 36;
[inflectionTrial2_1_36, meanProbabilities2_1_36, smoothedProbabilities2_1_36, upperBound2_1_36, lowerBound2_1_36] = analyzeLearningCurve(trials2(range2_1_36), responses2(range2_1_36, :));

% Analyze the second dataset for trials 36-end
range2_36_end = trials2 > 36;
[inflectionTrial2_36_end, meanProbabilities2_36_end, smoothedProbabilities2_36_end, upperBound2_36_end, lowerBound2_36_end] = analyzeLearningCurve(trials2(range2_36_end), responses2(range2_36_end, :));

% Plot the learning curves for trials 1-36
figure;

% Plot the first dataset (trials 1 to 36)
fill([trials1(range1_1_36); flipud(trials1(range1_1_36))], [upperBound1_1_36; flipud(lowerBound1_1_36)], [0.9 0.9 1], 'linestyle', 'none');
hold on;
plot(trials1(range1_1_36), smoothedProbabilities1_1_36, 'b-', 'LineWidth', 2);

% Mark the inflection point for the first dataset (trials 1 to 36)
if ~isnan(inflectionTrial1_1_36)
    plot(inflectionTrial1_1_36, smoothedProbabilities1_1_36(trials1(range1_1_36) == inflectionTrial1_1_36), 'ko', 'MarkerFaceColor', 'k');
    text(inflectionTrial1_1_36, smoothedProbabilities1_1_36(trials1(range1_1_36) == inflectionTrial1_1_36), ' Inflection', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end

% Plot the second dataset (trials 1 to 36)
fill([trials2(range2_1_36); flipud(trials2(range2_1_36))], [upperBound2_1_36; flipud(lowerBound2_1_36)], [1 0.9 0.9], 'linestyle', 'none');
plot(trials2(range2_1_36), smoothedProbabilities2_1_36, 'r-', 'LineWidth', 2);

% Mark the inflection point for the second dataset (trials 1 to 36)
if ~isnan(inflectionTrial2_1_36)
    plot(inflectionTrial2_1_36, smoothedProbabilities2_1_36(trials2(range2_1_36) == inflectionTrial2_1_36), 'ko', 'MarkerFaceColor', 'k');
    text(inflectionTrial2_1_36, smoothedProbabilities2_1_36(trials2(range2_1_36) == inflectionTrial2_1_36), ' Inflection', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end

xlabel('Trial Number (1 to 36)');
ylabel('Probability of Correct Response');
title('Learning Curves (Trials 1 to 36) with SEM Shading');
legend('First SEM', 'First Learning Curve', 'First Inflection', 'Second SEM', 'Second Learning Curve', 'Second Inflection');
grid on;
hold off;

% Plot the learning curves for trials 36-end
figure;

% Plot the first dataset (trials 36 to end)
fill([trials1(range1_36_end); flipud(trials1(range1_36_end))], [upperBound1_36_end; flipud(lowerBound1_36_end)], [0.9 0.9 1], 'linestyle', 'none');
hold on;
plot(trials1(range1_36_end), smoothedProbabilities1_36_end, 'b-', 'LineWidth', 2);

% Mark the inflection point for the first dataset (trials 36 to end)
if ~isnan(inflectionTrial1_36_end)
    plot(inflectionTrial1_36_end, smoothedProbabilities1_36_end(trials1(range1_36_end) == inflectionTrial1_36_end), 'ko', 'MarkerFaceColor', 'k');
    text(inflectionTrial1_36_end, smoothedProbabilities1_36_end(trials1(range1_36_end) == inflectionTrial1_36_end), ' Inflection', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end

% Plot the second dataset (trials 36 to end)
fill([trials2(range2_36_end); flipud(trials2(range2_36_end))], [upperBound2_36_end; flipud(lowerBound2_36_end)], [1 0.9 0.9], 'linestyle', 'none');
plot(trials2(range2_36_end), smoothedProbabilities2_36_end, 'r-', 'LineWidth', 2);

% Mark the inflection point for the second dataset (trials 36 to end)
if ~isnan(inflectionTrial2_36_end)
    plot(inflectionTrial2_36_end, smoothedProbabilities2_36_end(trials2(range2_36_end) == inflectionTrial2_36_end), 'ko', 'MarkerFaceColor', 'k');
    text(inflectionTrial2_36_end, smoothedProbabilities2_36_end(trials2(range2_36_end) == inflectionTrial2_36_end), ' Inflection', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end

xlabel('Trial Number (36 to End)');
ylabel('Probability of Correct Response');
title('Learning Curves (Trials 36 to End) with SEM Shading');
legend('First SEM', 'First Learning Curve', 'First Inflection', 'Second SEM', 'Second Learning Curve', 'Second Inflection');
grid on;
hold off;

% Display the inflection points
fprintf('Inflection point for the first dataset (trials 1 to 36): Trial %d\n', inflectionTrial1_1_36);
fprintf('Inflection point for the first dataset (trials 36 to end): Trial %d\n', inflectionTrial1_36_end);
fprintf('Inflection point for the second dataset (trials 1 to 36): Trial %d\n', inflectionTrial2_1_36);
fprintf('Inflection point for the second dataset (trials 36 to end): Trial %d\n', inflectionTrial2_36_end);
% Function to find the inflection point
function inflectionIndex = findInflectionPoint(data)
    secondDerivative = diff(diff(data));
    signChanges = diff(sign(secondDerivative));
    inflectionIndex = find(signChanges, 1);
end

% Function to analyze the learning curve
function [inflectionTrial, meanProbabilities, smoothedProbabilities, upperBound, lowerBound] = analyzeLearningCurve(trials, responses)
    numTrials = length(trials);
    numMice = size(responses, 2);
    predictedProbabilities = zeros(numTrials, numMice);

       for i = 1:numMice
        [b, ~, ~] = glmfit(trials, responses(:, i), 'binomial', 'link', 'logit');
        predictedProbabilities(:, i) = glmval(b, trials, 'logit');
       end
       
       
    meanProbabilities = mean(predictedProbabilities, 2);
    smoothedProbabilities = movmean(meanProbabilities, 5);
    SEM = std(predictedProbabilities, 0, 2) / sqrt(numMice);
    smoothedSEM = movmean(SEM, 5);

    upperBound = smoothedProbabilities + smoothedSEM;
    lowerBound = smoothedProbabilities - smoothedSEM;

    % Find the inflection point
    inflectionIndex = findInflectionPoint(smoothedProbabilities);
    if ~isempty(inflectionIndex)
        inflectionTrial = trials(inflectionIndex);
    else
        inflectionTrial = NaN;
    end
end


