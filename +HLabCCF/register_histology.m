function register_histology(data_path, image_folder)

processed_images_folder = fullfile(image_folder, 'processed');

% plane to view ('coronal', 'sagittal', 'transverse')
plane = 'coronal';

probe_save_name_suffix = '';

annotation_volume_location = fullfile(data_path, 'annotation_volume_10um_by_index.npy');
structure_tree_location = fullfile(data_path, 'structure_tree_safe_2017.csv');
template_volume_location = fullfile(data_path, 'template_volume_10um.npy');

%% GET PROBE TRAJECTORY POINTS

% load the reference brain and region annotations
if ~exist('av','var') || ~exist('st','var') || ~exist('tv','var')
    disp('loading reference atlas...')
    av = readNPY(annotation_volume_location);
    st = loadStructureTree(structure_tree_location);
    tv = readNPY(template_volume_location);
end

% select the plane for the viewer
if strcmp(plane,'coronal')
    av_plot = av;
    tv_plot = tv;
elseif strcmp(plane,'sagittal')
    av_plot = permute(av,[3 2 1]);
    tv_plot = permute(tv,[3 2 1]);
elseif strcmp(plane,'transverse')
    av_plot = permute(av,[2 3 1]);
    tv_plot = permute(tv,[2 3 1]);
end

%%
atlas_figure_browser = figure('Name', 'Atlas Viewer');
slice_figure_browser = figure('Name', 'Slice Viewer');

reference_size = size(tv_plot);
sliceBrowser(slice_figure_browser, processed_images_folder, atlas_figure_browser, reference_size);

% use application in Atlas Transform Viewer
% use this function if you have a processed_images_folder with appropriately processed .tif histology images
AtlasTransformBrowser(atlas_figure_browser, tv_plot, av_plot, st, slice_figure_browser, processed_images_folder, probe_save_name_suffix, plane);
end

