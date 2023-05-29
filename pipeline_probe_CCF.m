%% run this script section-by-section
clc
clear
addpath(genpath(fullfile(userpath, "npy-matlab", "npy-matlab")))
addpath(genpath(fullfile(userpath, "allenCCF")))
data_path = fullfile(userpath, "allenCCF", "data");

%% Step 1: load and adjust histology
image_folder = HLabCCF.load_histology();
pause; close all
HLabCCF.preprocess_histology(image_folder);
pause; close all

%% Step 2: register histology
% IMPORTANT: for probe points, the first one should be the very tip of the probe!!!
HLabCCF.register_histology(data_path, image_folder);
pause; close all

%% load preprocessed histology (dev)
image_folder = uigetdir(pwd,'Select folder with images for processing');

%% Step 3: process the clicked objects into CCF coordinates and save results
rez = HLabCCF.transform_clicked_points(image_folder);
for idx = 1:length(rez)
    item = rez{idx};
    save(fullfile(image_folder, 'processed', ['probe_points_transformed_probe_' int2str(idx) '.mat']), '-struct', 'item', '-v7.3');
end