function image_folder = load_histology(microns_per_pixel)
if nargin == 0
    microns_per_pixel = 1000 / 350; % 350px per mm
end

image_folder = uigetdir(pwd,'Select folder with images for processing');
save_folder = image_folder;

% name of images, in order anterior to posterior or vice versa
% once these are downsampled they will be named ['original name' '_processed.tif']
quest='Are the images .tif files?';
img=questdlg(quest,'Image');

if strcmpi(img, 'yes')
    image_file_names = dir([save_folder filesep '*.tif']); % get the contents of the image_folder
else
    image_type_prompt={'Enter the image file extension (ex: .tif, .jpg, etc.)'};
    type_title='Input';
    image_type=inputdlg(image_type_prompt,type_title);
    ext=strcat('*',image_type{1});
    image_file_names = dir([image_folder filesep ext]);
end

image_file_names = natsortfiles({image_file_names.name});

% if the images are individual slices (as opposed to images of multiple
% slices, which must be cropped using the cell CROP AND SAVE SLICES)
image_files_are_individual_slices = true;

% use images that are already at reference atlas resolution (here, 10um/pixel)
use_already_downsampled_image = microns_per_pixel == 10;

% pixel size parameters: microns_per_pixel of large images in the image
% folder (if use_already_downsampled_images is set to false);
% microns_per_pixel_after_downsampling should typically be set to 10 to match the atlas
% microns_per_pixel = 3.233;
microns_per_pixel_after_downsampling = 10;

% ----------------------
% additional parameters
% ----------------------

% increase gain if for some reason the images are not bright enough
gain = 1;


% finds or creates a folder location for processed images --
% a folder within save_folder called processed
folder_processed_images = fullfile(save_folder, 'processed');
if ~exist(folder_processed_images, 'dir')
    mkdir(folder_processed_images)
end

%%
histology_figure = figure('Name', 'Histology Viewer');
% Function to downsample and adjust histology image
HistologyBrowser(histology_figure, save_folder, image_folder, image_file_names, folder_processed_images, image_files_are_individual_slices, ...
    use_already_downsampled_image, microns_per_pixel, microns_per_pixel_after_downsampling, gain)

end

