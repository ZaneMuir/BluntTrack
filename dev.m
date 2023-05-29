%% 
addpath(genpath(fullfile(userpath, "npy-matlab", "npy-matlab")))
addpath(genpath(fullfile(userpath, "allenCCF")))
data_path = fullfile(userpath, "allenCCF", "data");


%% load and adjust histology
image_folder = HLabCCF.load_histology();
pause; close all
HLabCCF.preprocess_histology(image_folder);
pause; close all

%% register histology
HLabCCF.register_histology(data_path, image_folder);
pause; close all

%% process the clicked objects into CCF coordinates