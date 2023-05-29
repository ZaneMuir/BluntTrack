
%% Add to path

% addpath(genpath([pwd,'\Everyday_Functions']));
addpath(genpath([pwd,'/HistologicalReconstructionOfMultiERecordingSites/Kilosort_Related']));
addpath(genpath([pwd,'/HistologicalReconstructionOfMultiERecordingSites/npy-matlab-master']));
addpath(genpath([pwd,'/HistologicalReconstructionOfMultiERecordingSites/output']));
addpath(genpath([pwd,'/HistologicalReconstructionOfMultiERecordingSites/allenbrainatlascoords']));

try
    keep shp roi_locs
catch
end
close all
clc

dbstop if error

%% find associated files

% Kilosort files
kilosort_path=uigetdir(pwd,'Select folder containing kilosort files');
kilosort_files=dir(kilosort_path);
fpath = kilosort_path;

[~,kdate] = lastplstrtok(fpath,'\',1);
kdate = kdate(1:end-1);

% Metadata files
%metadata_path=uigetdir(pwd,'Select folder containing metadata');
metadata_path = [pwd,'\HistologicalReconstructionOfMultiERecordingSites\output\metadata\'];

% ROI files
%roi_path = uigetdir(pwd,'Select folder containing roi_location');
roi_path = [pwd,'\HistologicalReconstructionOfMultiERecordingSites\output\roi_location\'];

% Allen Atlas Locations
allen_path = [pwd,'\HistologicalReconstructionOfMultiERecordingSites\allenbrainatlascoords\'];

%% Set up roi_locations

% shp has Points, Alpha, HoleThreshold, and RegionThreshold
% The default alpha radius is a = criticalAlpha(shp,'all-points'), which is
%   the smallest alpha radius that produces an alpha shape enclosing all points. 
%   The extreme values of Alpha have the following conditions: If Alpha is Inf, 
%   then alphaShape produces the convex hull; If Alpha is 0, then the resulting 
%   alphaShape is empty. 
% Maximum area or volume of interior holes or voids to fill in, specified as a finite 
%   nonnegative scalar. For 3-D, HoleThreshold specifies the maximum volume of interior 
%   voids to fill in. Holes extending completely through the 3-D alpha shape cannot be filled in.
%   The default value is 0, so that alphaShape does not suppress any holes or voids. The application 
%   of the HoleThreshold and RegionThreshold properties is order-dependent. alphaShape fills in holes
%   before suppressing regions.
% Maximum area (2-D) or volume (3-D) of regions to suppress, specified as a finite nonnegative scalar.
%   The default value is 0, so that alphaShape does not suppress any regions. The application of the 
%   HoleThreshold and RegionThreshold properties is order-dependent. alphaShape fills in holes before 
%   suppressing regions


if(~exist('shp','var'))
    
    % dLGN
    load([allen_path,'dLGN\dlgn_apdvml_coord.mat']);
    dlgnshp = alphaShape(dlgn_apdvml(1,:)',dlgn_apdvml(2,:)',dlgn_apdvml(3,:)');
    %dlgnshp.Alpha = 0.0171; % 0.0166 from criticalAlpha, next value from alphaSpectrum
    
    % vLGN
    load([allen_path,'vLGN\vlgn_apdvml_coords.mat']);
    vlgnshp = alphaShape(apdvml_coord(1,:)',apdvml_coord(2,:)',apdvml_coord(3,:)');
    
    % IGL
    load([allen_path,'IGL\igl_apdvml_coords.mat']);
    iglshp = alphaShape(apdvml_coord(1,:)',apdvml_coord(2,:)',apdvml_coord(3,:)');
    
    % LP
    load([allen_path,'LP\lp_apdvml_coords.mat']);
    lpshp = alphaShape(apdvml_coord(1,:)',apdvml_coord(2,:)',apdvml_coord(3,:)');
    
    % MGd
    load([allen_path,'MGd\mgd_apdvml_coords.mat']);
    mgdshp = alphaShape(apdvml_coord(1,:)',apdvml_coord(2,:)',apdvml_coord(3,:)');
    
    % SGN
    load([allen_path,'SGn\sgn_apdvml_coords.mat']);
    sgnshp = alphaShape(apdvml_coord(1,:)',apdvml_coord(2,:)',apdvml_coord(3,:)');
    
    % TRN
    load([allen_path,'TRN\trn_apdvml_coords.mat']);
    trnshp = alphaShape(apdvml_coord(1,:)',apdvml_coord(2,:)',apdvml_coord(3,:)');
    
    % Hipp
    load([allen_path,'Hipp\hippf_apdvml_coords.mat']);
    hippshp = alphaShape(apdvml_coord(1,:)',apdvml_coord(2,:)',apdvml_coord(3,:)');
    
    % VPL
    load([allen_path,'VPL\vpl_apdvml_coords.mat']);
    vplshp = alphaShape(apdvml_coord(1,:)',apdvml_coord(2,:)',apdvml_coord(3,:)');
    
    % VPM
    load([allen_path,'VPM\vpm_apdvml_coords.mat']);
    vpmshp = alphaShape(apdvml_coord(1,:)',apdvml_coord(2,:)',apdvml_coord(3,:)');
    
    % PT - APN
    load([allen_path,'PT\APN\ptapn_apdvml_coords.mat']);
    ptapnshp = alphaShape(apdvml_coord(1,:)',apdvml_coord(2,:)',apdvml_coord(3,:)');
    
    % PT - MPA
    load([allen_path,'PT\MPA\ptmp_apdvml_coords.mat']);
    ptmpshp = alphaShape(apdvml_coord(1,:)',apdvml_coord(2,:)',apdvml_coord(3,:)');
    
    % PT - NOT
    load([allen_path,'PT\NOT\ptnot_apdvml_coords.mat']);
    ptnotshp = alphaShape(apdvml_coord(1,:)',apdvml_coord(2,:)',apdvml_coord(3,:)');
    
    % PT - NPC
    load([allen_path,'PT\NPC\ptnpc_apdvml_coords.mat']);
    ptnpcshp = alphaShape(apdvml_coord(1,:)',apdvml_coord(2,:)',apdvml_coord(3,:)');
    
    % PT - PPN
    load([allen_path,'PT\OPN\ptopn_apdvml_coords.mat']);
    ptopnshp = alphaShape(apdvml_coord(1,:)',apdvml_coord(2,:)',apdvml_coord(3,:)');
    
    % PT - PPN
    load([allen_path,'PT\PPN\ptppn_apdvml_coords.mat']);
    ptppnshp = alphaShape(apdvml_coord(1,:)',apdvml_coord(2,:)',apdvml_coord(3,:)');
    
    % SCs - optic
    load([allen_path,'SC\SCs\SCoptic\sco_apdvml_coords.mat']);
    scsoshp = alphaShape(apdvml_coord(1,:)',apdvml_coord(2,:)',apdvml_coord(3,:)');
    
    % SCs - superficial grey
    load([allen_path,'SC\SCs\SCsuperficialgrey\scsg_apdvml_coords.mat']);
    scssgshp = alphaShape(apdvml_coord(1,:)',apdvml_coord(2,:)',apdvml_coord(3,:)');
    
    % SCs - zonal
    load([allen_path,'SC\SCs\SCzonal\scz_apdvml_coords.mat']);
    scszshp = alphaShape(apdvml_coord(1,:)',apdvml_coord(2,:)',apdvml_coord(3,:)');
    
    % Combine it all together
    roi_locs = {'dLGN','vLGN','IGL','LP','MGd','SGn','TRN','Hipp','VPL','VPM','PTapn',...
        'PTmpa','PTnot','PTnpc','PTopn','PTppn','SCso','SCssg','SCsz'};
    shp = cell(1,length(roi_locs));
    
    shp{1} = dlgnshp;
    shp{2} = vlgnshp;
    shp{3} = iglshp;
    shp{4} = lpshp;
    shp{5} = mgdshp;
    shp{6} = sgnshp;
    shp{7} = trnshp;
    shp{8} = hippshp;
    shp{9} = vplshp;
    shp{10} = vpmshp;
    shp{11} = ptapnshp;
    shp{12} = ptmpshp;
    shp{13} = ptnotshp;
    shp{14} = ptnpcshp;
    shp{15} = ptopnshp;
    shp{16} = ptppnshp;
    shp{17} = scsoshp;
    shp{18} = scssgshp;
    shp{19} = scszshp;
end

%% Get the template depth

ksdir = loadKSdir(fpath,1);

% plotUnitScatter
clu = ksdir.clu(ismember(ksdir.clu,ksdir.cids));
spikeTemplates = ksdir.spikeTemplates(ismember(ksdir.clu,ksdir.cids));
thold = findTempForEachClu(clu, spikeTemplates);
tempPerClu = thold(~isnan(thold));

tempsUnW = zeros(size(ksdir.temps));
for t = 1:size(ksdir.temps,1)
    tempsUnW(t,:,:) = squeeze(ksdir.temps(t,:,:))*ksdir.winv;
end
tempChanAmps = squeeze(max(tempsUnW,[],2))-squeeze(min(tempsUnW,[],2));
tDepths = zeros(1,length(ksdir.cids));
cids = ksdir.cids;
for c = 1:length(ksdir.cids)
     thisTempInd = tempPerClu(c);
     thisTempAmps = tempChanAmps(thisTempInd+1,:);
     ampMed = median(thisTempAmps);
     peakYCtemp = ksdir.ycoords(thisTempAmps==max(thisTempAmps));
     tDepths(c) = peakYCtemp(1);
end

%% metadata collection

% select folder with the metadata
metadata_files=dir(metadata_path);
metadata_file_arr=string(1);
% loops through files to generate a list of the folder contents
for metadata_loop=1:size(metadata_files,1)
    metadata_name=strcat(metadata_path,'/',metadata_files(metadata_loop).name);
    if isfolder(metadata_name)==0
        metadata_file_arr=[metadata_file_arr;metadata_files(metadata_loop).name];
    end
end

km = contains(metadata_file_arr,kdate);
if(any(km)) % any of them equal to 1
    metadata_load = metadata_file_arr(km);
else
    [metadata_file_indx,~]=listdlg('ListString',metadata_file_arr);
    metadata_load=metadata_file_arr(metadata_file_indx);
end
load(metadata_load);  % load relevant file

[a,~]=size(metadata_struct);
penetrationNumbers=a;

%% Load in ROI file
% roi_files = dir(roi_path);
% roi_file_arr=string(1);
% % loops through files to generate a list of the folder contents
% for roi_loop=1:size(roi_files,1)
%     roi_name=strcat(roi_path,'/',roi_files(roi_loop).name);
%     if isfolder(roi_name)==0
%         roi_file_arr=[roi_file_arr;roi_files(roi_loop).name];
%     end
% end
% 
% kr =  contains(roi_file_arr,kdate);
% if(any(kr)) % any of them equal to 1
%     roi_load = roi_file_arr(kr);
% else
%     [roi_file_indx,~]=listdlg('ListString',roi_file_arr);
%     roi_load=roi_file_arr(roi_file_indx);
% end
% load(roi_load);  % load relevant file

%% Begin to go through the penetrations

ks_metadata_struct(a).cluster_ids = [];

for penetrations=1:a
      
    % user input for bottom most / top most conductor
    quest='Is 0 the top or bottom most conductor?';
    conductor_position=questdlg(quest,'Conductors','Top','Bottom','Top');
    
    %% find cluster coordinates
    t=struct2table(metadata_struct);
    conductors=t.(5);
    
    if strcmpi(conductor_position,'Top')
        zPoint=size(conductors,1);
        lPoint = 1;
    elseif strcmpi(conductor_position,'Bottom')
        zPoint=1;
        lPoint = size(conductors,1);
    end
    
    topCond=conductors(zPoint,:);
    tDepthsMm=tDepths/1000;
    
    % determining angles in radians use trig
    track_index=penetrations;
    %roi_struct_correct = roi_struct;
    roi_struct = metadata_struct; % she forgot the damn conversion; also metadata is simialr to roi_struct just with additional stuff
    
%     [b,~]=size(roi_struct(track_index).zfit);
%     if b==1
%         [~,b]=size(roi_struct(track_index).zfit);
%     end
%     zdiff=(roi_struct(track_index).zfit(b)-roi_struct(track_index).zfit(1)); % find distance between points to calculate angles
%     xdiff=(roi_struct(track_index).xfit(b)-roi_struct(track_index).xfit(1));
%     ydiff=(roi_struct(track_index).yfit(b)-roi_struct(track_index).yfit(1));
    zdiff = conductors(lPoint,3)-conductors(zPoint,3);
    xdiff = conductors(lPoint,1)-conductors(zPoint,1);
    ydiff = conductors(lPoint,2)-conductors(zPoint,2);
    
    x_z_rad=atan(zdiff/xdiff);
    r=sqrt((zdiff^2)+(xdiff^2));
    y_rad=atan(ydiff/r);
    
    clus_coords=zeros(length(tDepthsMm),3);
    for clus_ind=1:size(tDepthsMm,2)
        
        y_adjust=tDepthsMm(clus_ind)*sin(y_rad);
        r_adjust=y_adjust/tan(y_rad);
        x_adjust=r_adjust*cos(x_z_rad);
        z_adjust=r_adjust*sin(x_z_rad);
        
        clus_x=topCond(1)+x_adjust;
        clus_y=topCond(2)+y_adjust;
        clus_z=topCond(3)+z_adjust;
        
        clus_coords(clus_ind,:)=[clus_x,clus_y,clus_z];   
    end
    
    % Determine placement of the clusters (reorganize by tDepths)
    %[~,ti] = sort(tDepths);
    
    clus_coords_temp = [clus_coords(:,1) clus_coords(:,3) clus_coords(:,2)];
    %clus_coords_temp = clus_coords_temp(ti,:);
    
    areacheck = zeros(size(clus_coords,1),length(roi_locs));
    nncheck = zeros(size(clus_coords,1),length(roi_locs),3);
    nneuc = zeros(size(clus_coords,1),length(roi_locs));
    
    for s = find(~cellfun('isempty',shp))
        currshp = shp{s};
        
		%tri = alphaTriangulation(currshp);
		%areacheck(:,s) = pointLocation(tri,clus_coords_temp);
		
        areacheck(:,s) = inShape(currshp,clus_coords_temp);
        
        [nni,nneuc(:,s)] = nearestNeighbor(currshp,clus_coords_temp);
        nncheck(:,s,:) = currshp.Points(nni,:);
    end
    
    nneuc(nneuc==0) = inf;
    [~,minneuci] = min(nneuc,[],2);
    
    % Just look at all together in a table
    % cluster number, cluster coordinates, area check for each area
        
    areacolumn = cell(size(areacheck,1),1);
    areacolumn(cellfun('isempty',areacolumn)) = {NaN};
    for ar = 1:size(areacheck,2)
        currlog = areacheck(:,ar);
        areacolumn(logical(currlog)) = roi_locs(ar);
    end
    
    [row,col] = find(areacheck==1);
    
    %nneuccolumn = cell(length(minneuci),1);
    %nneuccolumn(cellfun('isempty',nneuccolumn)) = {NaN};
    %ac = sum(areacheck,2); 
    ntemp = nneuc; 
    for acc = 1:length(minneuci)
        %if(ac(acc)<=0)
        if(~ismember(acc,row))
            continue;
        else
            ntemp(acc,col) = NaN;
        end
    end
    [minneuc,minneucit] = min(ntemp,[],2);
    for n = 1:length(minneucit)
        nneuccolumn(n) = roi_locs(minneucit(n));
    end
    
    T = table(cids',tDepths',clus_coords_temp,areacolumn,nneuccolumn',minneuc,...
        'VariableNames',{'ClusterID','ElecDepth','ClusterCoordsApDvMl','ClusterROI','NearestROItoClusterROI','EucDistinMM'});
       
end

disp(T)

%% Additional functions

function spikeStruct = loadKSdir(ksDir, excludeNoise)

% if ~isempty(varargin)
%     params = varargin{1};
% else
     params = [];
% end

params.excludeNoise = excludeNoise;

% if ~isfield(params, 'excludeNoise')
%     params.excludeNoise = true;
% end
% if ~isfield(params, 'loadPCs')
params.loadPCs = false;
% end

% load spike data

spikeStruct = loadParamsPy(fullfile(ksDir, 'params.py'));

ss = readNPY(fullfile(ksDir, 'spike_times.npy'));
st = double(ss)/spikeStruct.sample_rate;
spikeTemplates = readNPY(fullfile(ksDir, 'spike_templates.npy')); % note: zero-indexed

if exist(fullfile(ksDir, 'spike_clusters.npy'),'file')
    clu = readNPY(fullfile(ksDir, 'spike_clusters.npy'));
else
    clu = spikeTemplates;
end

tempScalingAmps = readNPY(fullfile(ksDir, 'amplitudes.npy'));

if params.loadPCs
    pcFeat = readNPY(fullfile(ksDir,'pc_features.npy')); % nSpikes x nFeatures x nLocalChannels
    pcFeatInd = readNPY(fullfile(ksDir,'pc_feature_ind.npy')); % nTemplates x nLocalChannels
else
    pcFeat = [];
    pcFeatInd = [];
end

cgsFile = '';
if exist(fullfile(ksDir, 'cluster_groups.csv'),'file') 
    cgsFile = fullfile(ksDir, 'cluster_groups.csv');
end
if exist(fullfile(ksDir, 'cluster_group.tsv'),'file') 
   cgsFile = fullfile(ksDir, 'cluster_group.tsv');
end 
if ~isempty(cgsFile)
    [cids, cgs] = readClusterGroupsCSV(cgsFile);

    if params.excludeNoise
        noiseClusters = cids(cgs==0);

        st = st(~ismember(clu, noiseClusters));
        spikeTemplates = spikeTemplates(~ismember(clu, noiseClusters));
        tempScalingAmps = tempScalingAmps(~ismember(clu, noiseClusters));        
        
        if params.loadPCs
            pcFeat = pcFeat(~ismember(clu, noiseClusters), :,:);
            %pcFeatInd = pcFeatInd(~ismember(cids, noiseClusters),:);
        end
        
        clu = clu(~ismember(clu, noiseClusters));
        cgs = cgs(~ismember(cids, noiseClusters));
        cids = cids(~ismember(cids, noiseClusters));
        
        
    end
    
else
    clu = spikeTemplates;
    
    cids = unique(spikeTemplates);
    cgs = 3*ones(size(cids));
end
    

coords = readNPY(fullfile(ksDir, 'channel_positions.npy'));
ycoords = coords(:,2); xcoords = coords(:,1);
temps = readNPY(fullfile(ksDir, 'templates.npy'));

winv = readNPY(fullfile(ksDir, 'whitening_mat_inv.npy'));

spikeStruct.st = st;
spikeStruct.spikeTemplates = spikeTemplates;
spikeStruct.clu = clu;
spikeStruct.tempScalingAmps = tempScalingAmps;
spikeStruct.cgs = cgs;
spikeStruct.cids = cids;
spikeStruct.xcoords = xcoords;
spikeStruct.ycoords = ycoords;
spikeStruct.temps = temps;
spikeStruct.winv = winv;
spikeStruct.pcFeat = pcFeat;
spikeStruct.pcFeatInd = pcFeatInd;

end

function [b4str,aftstr] = lastplstrtok(YourString,delimiter,n)

if(nargin~=3)
    error('Need 2 inputs');
end

if(~ischar(YourString))
    error('Please enter a string to search');
end

if(~ischar(delimiter))
    error('Please enter a string to delimit');
end

lastdot_pos = find(YourString == delimiter, n, 'last');
lastdot_pos = lastdot_pos(1);
b4str = YourString(1:lastdot_pos);
aftstr = YourString(lastdot_pos+1:end);

end