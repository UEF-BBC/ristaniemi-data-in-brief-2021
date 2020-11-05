% Example script for loading NIRS ligament data.
%
% The script visualizes the average spectrum of each ligament type. Make sure
% the file "Data_NIRS_measurements.xlsx" is in the same directory before running
% the script.
%
% Jari Torniainen
% Department of Applied Physics, Univeristy of Eastern Finland, 2019

clear all; close all; clc;

% --- Load data ---
data_vis = readtable('Data_NIRS_measurements.xlsx', 'readrownames', true, 'sheet', 'AvaSpec-ULS2048L-USB2');
data_nir = readtable('Data_NIRS_measurements.xlsx', 'readrownames', true, 'sheet', 'AvaSpec-NIR256-2.5-HSC');

% NOTE: You can read EXACT wavelength values from the column names of the data. I am taking a shortcut by
% generating them with the linspace-function.
wave_vis = linspace(319.43, 1093.04, 1343);
wave_nir = linspace(943.84, 2497.4, 243);

% --- Plot ligament type averages ---
ligament_types = {'ACL', 'PCL', 'LCL', 'MCL', 'PT'};
figure();
for ligament_type = ligament_types
    mask = contains(data_vis.Properties.RowNames, ligament_type);
    mean_spectrum_vis = mean(data_vis{mask, 2:end}, 1);
    subplot(1, 2, 1)
    hold on;
    plot(wave_vis, mean_spectrum_vis)

    mask = contains(data_nir.Properties.RowNames, ligament_type);
    mean_spectrum_nir = mean(data_nir{mask, 2:end}, 1);
    subplot(1, 2, 2)
    hold on;
    plot(wave_nir, mean_spectrum_nir)
end

subplot(1, 2, 1)
legend(ligament_types);
ylabel('Absorbance');
xlabel('Wavelength [nm]');
title('AvaSpec-ULS2048L-USB2');

subplot(1, 2, 2)
xlabel('Wavelength [nm]');
title('AvaSpec-NIR256-2.5-HSC');
