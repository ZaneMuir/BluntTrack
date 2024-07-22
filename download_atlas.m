% http://data.cortexlab.net/allenCCF/annotation_volume_10um_by_index.npy
clear
clc
%%
data_path = fullfile(userpath, "allenCCF");
if ~exist(data_path, "dir")
    mkdir(data_path);
end

%% annotation_volume_10um
filename = fullfile(data_path, 'annotation_volume_10um_by_index.npy');
if ~exist(filename, "file")
    disp("dowloading annotation atlas...")
    % websave(filename, "http://data.cortexlab.net/allenCCF/annotation_volume_10um_by_index.npy");
    websave(filename, "https://figshare.com/ndownloader/files/44925493")
    disp("annotation atlas downloaded.")
else
    disp("annotation atlas exists.")
end

%% template_volume_10um
filename = fullfile(data_path, 'template_volume_10um.npy');
if ~exist(filename, "file")
    disp("dowloading template atlas...")
    % websave(filename, "http://data.cortexlab.net/allenCCF/template_volume_10um.npy");
    websave(filename, "https://figshare.com/ndownloader/files/44925496")
    disp("template atlas downloaded.")
else
    disp("template atlas exists.")
end

%% structure_tree_safe_2017
filename = fullfile(data_path, 'structure_tree_safe_2017.csv');
if ~exist(filename, "file")
    % disp("dowloading structure tree 2017...")
    % websave(filename, "http://data.cortexlab.net/allenCCF/structure_tree_safe_2017.csv");
    % disp("structure tree 2017 downloaded.")
    disp("structure_tree_safe_2017.csv is missing, try to sync with github.")
else
    disp("structure tree 2017 exists.")
end

%% structure_tree_safe
filename = fullfile(data_path, 'structure_tree_safe.csv');
if ~exist(filename, "file")
    % disp("dowloading structure tree...")
    % websave(filename, "http://data.cortexlab.net/allenCCF/structure_tree_safe.csv");
    % disp("structure tree downloaded.")
    disp("structure_tree_safe.csv is missing, try to sync with github.")
else
    disp("structure tree exists.")
end

%%
