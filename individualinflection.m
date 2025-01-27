clc;
close all;

% --- Load the first dataset ---
filename1 = 'C:\Users\Anacker1\Desktop\logistic regression\vCA1-ACQ-INH-Cont.xlsx'; % Update with your actual file path
opts1 = detectImportOptions(filename1);
opts1.DataRange = 'A2';  % Start importing data from the second row
data1 = readmatrix(filename1, opts1);
trials1 = data1(:, 1);
responses1 = data1(:, 2:end);

% --- Load the second dataset ---
filename2 = 'C:\Users\Anacker1\Desktop\logistic regression\vCA1-ACQ-INH-Dreadd.xlsx'; % Update with your actual file path
opts2 = detectImportOptions(filename2);
opts2.DataRange = 'A2';  % Start importing data from the second row
data2 = readmatrix(filename2, opts2);
trials2 = data2(:, 1);
responses2 = data2(:, 2:end);

% --- Plot inflection points for each Control mouse in two ranges ---
fprintf('Inflection points for each Control mouse:\n');
figure;
numMice1 = size(responses1, 2); % Number of mice in the first dataset

for i = 1:numMice1
    % --- Analyze the first range (trials 1–36) ---
    range1_1_36 = trials1 >= 1 & trials1 <= 36;
    [inflectionTrial1_1_36, meanProbabilities1_1_36, smoothedProbabilities1_1_36, upperBound1_1_36, lowerBound1_1_36] = analyzeLearningCurve(trials1(range1_1_36), responses1(range1_1_36, i));

    % Plot the regression curve for the first range
    subplot(numMice1, 2, 2*i-1);
    fill([trials1(range1_1_36); flipud(trials1(range1_1_36))], [upperBound1_1_36; flipud(lowerBound1_1_36)], [0.9 0.9 1], 'linestyle', 'none');
    hold on;
    plot(trials1(range1_1_36), smoothedProbabilities1_1_36, 'b-', 'LineWidth', 2);
    
    % Mark the inflection point for trials 1–36
    if ~isnan(inflectionTrial1_1_36)
        plot(inflectionTrial1_1_36, smoothedProbabilities1_1_36(trials1(range1_1_36) == inflectionTrial1_1_36), 'ko', 'MarkerFaceColor', 'k');
        text(inflectionTrial1_1_36, smoothedProbabilities1_1_36(trials1(range1_1_36) == inflectionTrial1_1_36), ' Inflection', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    end
    
    xlabel('Trial Number (1 to 36)');
    ylabel('Probability of Correct Response');
    title(sprintf('Mouse %d (Control, 1–36)', i));
    grid on;
    hold off;
    
    % --- Analyze the second range (trials 36–end) ---
    range1_36_end = trials1 > 36;
    [inflectionTrial1_36_end, meanProbabilities1_36_end, smoothedProbabilities1_36_end, upperBound1_36_end, lowerBound1_36_end] = analyzeLearningCurve(trials1(range1_36_end), responses1(range1_36_end, i));
    
    % Plot the regression curve for the second range
    subplot(numMice1, 2, 2*i);
    fill([trials1(range1_36_end); flipud(trials1(range1_36_end))], [upperBound1_36_end; flipud(lowerBound1_36_end)], [0.9 0.9 1], 'linestyle', 'none');
    hold on;
    plot(trials1(range1_36_end), smoothedProbabilities1_36_end, 'b-', 'LineWidth', 2);
    
    % Mark the inflection point for trials 36–end
    if ~isnan(inflectionTrial1_36_end)
        plot(inflectionTrial1_36_end, smoothedProbabilities1_36_end(trials1(range1_36_end) == inflectionTrial1_36_end), 'ko', 'MarkerFaceColor', 'k');
        text(inflectionTrial1_36_end, smoothedProbabilities1_36_end(trials1(range1_36_end) == inflectionTrial1_36_end), ' Inflection', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    end
    
    xlabel('Trial Number (36 to end)');
    ylabel('Probability of Correct Response');
    title(sprintf('Mouse %d (Control, 36–end)', i));
    grid on;
    hold off;
    
    % Print inflection points for both ranges for each mouse
    fprintf('Mouse %d (Control): Inflection at trial %d (1-36), and trial %d (36-end)\n', i, inflectionTrial1_1_36, inflectionTrial1_36_end);
end

% --- Plot inflection points for each DREADD mouse in two ranges ---
fprintf('\nInflection points for each DREADD mouse:\n');
figure;
numMice2 = size(responses2, 2); % Number of mice in the second dataset

for i = 1:numMice2
    % --- Analyze the first range (trials 1–36) ---
    range2_1_36 = trials2 >= 1 & trials2 <= 36;
    [inflectionTrial2_1_36, meanProbabilities2_1_36, smoothedProbabilities2_1_36, upperBound2_1_36, lowerBound2_1_36] = analyzeLearningCurve(trials2(range2_1_36), responses2(range2_1_36, i));

    % Plot the regression curve for the first range
    subplot(numMice2, 2, 2*i-1);
    fill([trials2(range2_1_36); flipud(trials2(range2_1_36))], [upperBound2_1_36; flipud(lowerBound2_1_36)], [1 0.9 0.9], 'linestyle', 'none');
    hold on;
    plot(trials2(range2_1_36), smoothedProbabilities2_1_36, 'r-', 'LineWidth', 2);
    
    % Mark the inflection point for trials 1–36
    if ~isnan(inflectionTrial2_1_36)
        plot(inflectionTrial2_1_36, smoothedProbabilities2_1_36(trials2(range2_1_36) == inflectionTrial2_1_36), 'ko', 'MarkerFaceColor', 'k');
        text(inflectionTrial2_1_36, smoothedProbabilities2_1_36(trials2(range2_1_36) == inflectionTrial2_1_36), ' Inflection', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    end
    
    xlabel('Trial Number (1 to 36)');
    ylabel('Probability of Correct Response');
    title(sprintf('Mouse %d (DREADD, 1–36)', i));
    grid on;
    hold off;
    
    % --- Analyze the second range (trials 36–end) ---
    range2_36_end = trials2 > 36;
    [inflectionTrial2_36_end, meanProbabilities2_36_end, smoothedProbabilities2_36_end, upperBound2_36_end, lowerBound2_36_end] = analyzeLearningCurve(trials2(range2_36_end), responses2(range2_36_end, i));
    
    % Plot the regression curve for the second range
    subplot(numMice2, 2, 2*i);
    fill([trials2(range2_36_end); flipud(trials2(range2_36_end))], [upperBound2_36_end; flipud(lowerBound2_36_end)], [1 0.9 0.9], 'linestyle', 'none');
    hold on;
    plot(trials2(range2_36_end), smoothedProbabilities2_36_end, 'r-', 'LineWidth', 2);
    
       % Mark the inflection point for trials 36–end
    if ~isnan(inflectionTrial2_36_end)
        plot(inflectionTrial2_36_end, smoothedProbabilities2_36_end(trials2(range2_36_end) == inflectionTrial2_36_end), 'ko', 'MarkerFaceColor', 'k');
        text(inflectionTrial2_36_end, smoothedProbabilities2_36_end(trials2(range2_36_end) == inflectionTrial2_36_end), ' Inflection', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    end
    
    xlabel('Trial Number (36 to end)');
    ylabel('Probability of Correct Response');
    title(sprintf('Mouse %d (DREADD, 36–end)', i));
    grid on;
    hold off;
    xlabel('Trial Number (36 to end)')
    
    % Print inflection points for both ranges for each DREADD mouse
    fprintf('Mouse %d (DREADD): Inflection at trial %d (1-36), and trial %d (36-end)\n', i, inflectionTrial2_1_36, inflectionTrial2_36_end);
end

% --- Function Definitions ---

% Function to find the inflection point
function inflectionIndex = findInflectionPoint(data)
    secondDerivative = diff(diff(data));
    signChanges = diff(sign(secondDerivative));
    inflectionIndex = find(signChanges, 1);
end

% Function to analyze the learning curve
function [inflectionTrial, meanProbabilities, smoothedProbabilities, upperBound, lowerBound] = analyzeLearningCurve(trials, responses)
    numTrials = length(trials);
    [numDataPoints, numMice] = size(responses);  % For each individual mouse's data
    
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
