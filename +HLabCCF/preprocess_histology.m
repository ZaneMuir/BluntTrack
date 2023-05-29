function preprocess_histology(image_folder)

    folder_processed_images = fullfile(image_folder, 'processed');
    
    % plane to view ('coronal', 'sagittal', 'transverse')
    plane = 'coronal';
    
    % size in pixels of reference atlas brain. For coronal slice, this is 800 x 1140
    if strcmp(plane,'coronal')
        atlas_reference_size = [800 1140]; 
    elseif strcmp(plane,'sagittal')
        atlas_reference_size = [800 1320]; 
    elseif strcmp(plane,'transverse')
        atlas_reference_size = [1140 1320];
    end

    slice_figure = figure('Name','Slice Viewer');
    SliceFlipper(slice_figure, folder_processed_images, atlas_reference_size)
end

