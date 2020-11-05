% Example script for constructing regression models with NIRS and reference 
% variables.
%
% This script attempts to predict the tissue water content from the NIR
% spectrum with partial least squares regression. Model performance is
% estimated with 9 fold cross-validation. The used modeling approach
% was kept simple to favour readability of the example.
%
% Jari Torniainen
% Department of Applied Physics, Univeristy of Eastern Finland, 2019

rng(6700);

% Load reference variable (proportional water content)
opts = detectImportOptions('Data_biomechanics_and_biochemistry.xlsx')
opts.VariableNamesRange = 'A11:AE11';
opts.DataRange = 'A13:AE62';
data_ref = readtable('Data_biomechanics_and_biochemistry.xlsx', opts);
water_content = data_ref{:, 25};

% Load NIR spectra
data_nir = readtable('Data_NIRS_measurements.xlsx', 'readrownames', true, 'sheet', 'AvaSpec-NIR256-2.5-HSC');
wave_nir = linspace(943.84, 2497.4, 243);

% Matching the spectra to reference variables. There are five NIRS measurements per each reference 
% variable so we just average them here.
ligament_types = {'ACL', 'PCL', 'MCL', 'LCL', 'PT'};
nirs = [];

for m = 1:10
    for ligament_type = ligament_types
        pattern = sprintf('M%02d_%s', m, ligament_type{1});
        mask = contains(data_nir.Properties.RowNames, pattern);
        spectrum = mean(data_nir{mask, 2:end}, 1);
        nirs(end + 1, :) = spectrum;
    end
end

% Limiting the full spectrum to a subset around the water peak of 1400 - 1500 nm
wave_min = 1000;
wave_max = 1920;
mask = wave_nir >= wave_min & wave_nir <= wave_max;
wave_nir = wave_nir(mask);
nirs = nirs(:, mask);

nirs = nirs'; % Each spectrum is flipped to a column vector

nirs = diff(nirs); % First order spectral derivative
nirs = sgolayfilt(nirs, 2, 61); % Savitzky-Golay filtering (3rd order polynomial, 61 points)

wave_nir = linspace(wave_min, wave_max, size(nirs, 1));

% Plot the preprocessed set of spectra
figure()
plot(wave_nir, nirs)
xlabel('Wavelength [nm]');
ylabel('Absorbance');
title('Filtered 1st derivative');

% Initialize 9 fold cross-validation
n_folds = 9;
cv = cvpartition(50, 'kfold', n_folds);

% For each fold: construct and test a PLRS model
figure();
for fold = 1:n_folds
    i_train = cv.training(fold);
    i_test = cv.test(fold);

    X_train = nirs(:, i_train)';
    X_test = nirs(:, i_test)';

    y_train = water_content(i_train);
    y_test = water_content(i_test);

    [~, ~, ~, ~, beta_values, pct_var, msep, stats] = plsregress(X_train, y_train, 3);

    y_pred = [ones(size(X_test, 1), 1) X_test] * beta_values;

    % Calculate correlation and root mean squared error as performance metrics
    r = corrcoef(y_pred, y_test);
    r = r(1, 2);
    rmse = sqrt(mean((y_test - y_pred).^2));

    % Plot expected predicted values as a scatter plot. Also print performance metrics
    subplot(3, 3, fold);
    plot(y_pred, y_test, 'ko');
    axis('square');
    xlim([60, 90]);
    ylim([60, 90]);
    xticks([60, 65, 70, 75, 80, 85, 90]);
    title(sprintf('fold=%d, r=%0.2f, RMSE=%0.2f%%',fold, r, rmse))
    xlabel('Predicted');
    ylabel('True');
    grid();
end
