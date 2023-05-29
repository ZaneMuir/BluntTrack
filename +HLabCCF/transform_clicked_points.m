function roi_location = transform_clicked_points(image_folder)

%% set parameters and load data
save_folder = fullfile(image_folder, 'processed');

% annotation_volume_location = fullfile(data_path, 'annotation_volume_10um_by_index.npy');
% structure_tree_location = fullfile(data_path, 'structure_tree_safe_2017.csv');
% template_volume_location = fullfile(data_path, 'template_volume_10um.npy');

% name of the saved object points
object_save_name_suffix = '';

% either set to 'all' or a list of indices from the clicked objects in this file, e.g. [2,3]
objects_to_analyze = 'all';

% plane used to view when points were clicked ('coronal' -- most common, 'sagittal', 'transverse')
plane = 'coronal';

% brain figure black or white
% black_brain = true;

% load the reference brain annotations
%     if ~exist('av','var') || ~exist('st','var')
%         disp('loading reference atlas...')
%         av = readNPY(annotation_volume_location);
%         st = loadStructureTree(structure_tree_location);
%     end


% load object points
objectPoints = load(fullfile(save_folder, ['probe_points' object_save_name_suffix]));

% determine which objects to analyze
if strcmp(objects_to_analyze,'all')
    objects = 1:size(objectPoints.pointList.pointList,1);
else
    objects = objects_to_analyze;
end

%% BRING UP THE RELEVANT DATA FOR EACH PROBE POINTS, FOR FURTHER ANALYSIS

% initialize cell array containing info on each clicked point
if length(objects) > 1
    %roi_annotation = cell(length(objects),1);
    roi_location = cell(length(objects),1);
end

% generate needed values
bregma = allenCCFbregma(); % bregma position in reference data space
atlas_resolution = 0.010; % mm

for object_num = objects
    selected_object = objects(object_num);

    % get the object points for the currently analyzed object
    if strcmp(plane,'coronal')
        curr_objectPoints = objectPoints.pointList.pointList{selected_object,1}(:, [3 2 1]);
    elseif strcmp(plane,'sagittal')
        curr_objectPoints = objectPoints.pointList.pointList{selected_object,1}(:, [1 2 3]);
    elseif strcmp(plane,'transverse')
        curr_objectPoints = objectPoints.pointList.pointList{selected_object,1}(:, [1 3 2]);
    end

    % use the point's position in the atlas to get the AP, DV, and ML coordinates
    ap = -(curr_objectPoints(:,1)-bregma(1))*atlas_resolution;
    dv = (curr_objectPoints(:,2)-bregma(2))*atlas_resolution;
    ml = (curr_objectPoints(:,3)-bregma(3))*atlas_resolution;

    roi_location_curr = [ap dv ml];
    probe_tip = roi_location_curr(1, :);
    probe_vec = mean(diff(roi_location_curr, 1, 1), 1);
    probe_uvec = probe_vec / norm(probe_vec);
    roi_location{object_num} = struct('ROI', roi_location_curr, ...
        'probe_tip', probe_tip, 'probe_uvec', probe_uvec, ...
        'cooridnates', ["anterior", "ventral", "right"], ...
        'unit', "mm");
end
end

