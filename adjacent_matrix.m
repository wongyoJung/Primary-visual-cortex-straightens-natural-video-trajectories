close all;
clear;

addpath(genpath('./library'));

% conditions
NUM_SCALE               = 2;
NUM_CATEGORIES          = 3;
NUM_MOVIES              = 10;
NUM_FRAMES              = 11;
NUM_POPULATIONS         = 9;

STIM_IMAGE_ON_SEC       = 0.2;
STIM_IMAGE_OFF_SEC      = 0.1;

% enumerations
NATURAL_TYPE_ENUM       = 1;
SYNTHETIC_TYPE_ENUM     = 2;
CONTRAST_TYPE_ENUM      = 3; % N/A


% LOAD DATA
% 1. results from curvature estimation
load(fullfile('./curvatureEstimations/', 'results_combined_table_235x.mat'), 'Acutedatacombined235x');

% 2. results from image computable model
agg_filepath            = fullfile('modelResults/', 'aggregated_LNLN_nat_sigG2_sqrt.mat');
load(agg_filepath, ...
    'c_V1_pixel', 'c_V1_model', ...
    'd_V1_pixel', 'd_V1_model', ...
    'corr_V1_model_data', ...
    'monkey_id_list', ...
    'penetration_id_list', 'num_units_recording', ...
    'natural_movie_labels', 'artificial_movie_labels', ...
    'weights_sd');
delta_curvature         = rad2deg(c_V1_model - c_V1_pixel);
distance_model          = d_V1_model;

% 3. data quality (in terms of how predictable one half is to the other)
DATA_FOLDER     = 'V1_responses/';
load(fullfile(DATA_FOLDER, 'data_correlations.mat'), 'data_correlations');

% 4. load recordings from image sequences (spike counts)
recordings      = dir( fullfile(DATA_FOLDER, '*_image_sequences.mat'));
frames_dim      = 4;
repeats_dim     = 6;

% 5. load bootstraps
load(fullfile(DATA_FOLDER, 'bootstrap_table.mat'), 'bootstrap_table');

categories              = [NATURAL_TYPE_ENUM, SYNTHETIC_TYPE_ENUM];

% list of categories
category_list           = Acutedatacombined235x.sequencetype;

% list of dither (step sizes)
dither_list             = Acutedatacombined235x.dither;

% list of recordings
recordingLabel          = Acutedatacombined235x.Dataset;
uniq_label_cat          = unique(Acutedatacombined235x.Dataset);
num_recordings          = numel(uniq_label_cat);
uniq_labels             = cell(numel(uniq_label_cat, 1));

typeofzooms = 2;
categories = 3;
movies = 10;
frames = 11;

%% adjacent matrix - cross correlation between pairs of neurons
angles = nan(num_recordings,typeofzooms,categories,movies,frames);
for i = 1:num_recordings
    
    % load spike data for each table entry
    i_recording_label           = strtok(recordings(i).name, '.');
    i_recording_label_mask      = i_recording_label == recordingLabel;
    i_filepath                  = fullfile(recordings(i).folder, recordings(i).name);
    load(i_filepath, 'naturalStruct');
    
    neuralresponse = naturalStruct.sortedResponses;    

    for z=1:typeofzooms
    for c=1:2
    for m=1:movies

            nanrows = [];
            % remove repeats having nan
            for f=1:frames
                cell_resps = squeeze(neuralresponse(z,c,m,f,:,:));
                [nanrow,nancolumn] = find(isnan(cell_resps));
                
                for k=1:length(nanrow)
                    row = nanrow(k);
                    col = nancolumn(k);
                    
                    meanofunit = mean(neuralresponse(z,c,m,f,row,:),'omitnan');
                    neuralresponse(z,c,m,f,row,col) = meanofunit;
                end
            end
            cell_res_allframe = neuralresponse(z,c,m,:,:,:);
            cell_res_allframe = squeeze(cell_res_allframe);
            % get correlation coeffecient
            numofunits = size(cell_res_allframe,2);
            coeffs = zeros(numofunits);
            for f=1:frames
                cell_res = cell_res_allframe(f,:,:);
                cell_res = squeeze(cell_res);
                coefff = corrcoef(cell_res');
                coeffs = coeffs+coefff;
            end
            coeffs = coeffs./frames;
            figure; imagesc(coeffs); colorbar;

%             for f=1:frames-1
%                 cell_res = cell_res_allframe(f,:,:);
%                 cell_res = squeeze(cell_res);
%                 std_cell_res = std(cell_res,[],2);
% 
%                 next_cell_res = cell_res_allframe(f+1,:,:);
%                 next_cell_res = squeeze(next_cell_res);
%                 std_cell_res_next = std(next_cell_res,[],2);
%                 
%                 coeff_val = corrcoef(std_cell_res,std_cell_res_next);
%                 angles(i,z,c,m,f) = coeff_val(1,2);
%             end

    end
    end
    end
end


% natural images
natural_frames_mean_arr = [];
natural_frames_sem_arr = [];
for f=1:frames
    
    reshaped_angles = reshape(angles(:,:,1,:,f),[],1);
    
    natural_frames_mean = mean(reshaped_angles,'omitnan');
    natural_frames_std = std(reshaped_angles,'omitnan');
    natural_frames_sem = natural_frames_std / sqrt(length(reshaped_angles));
    
    natural_frames_mean_arr = [natural_frames_mean_arr; natural_frames_mean];
    natural_frames_sem_arr = [natural_frames_sem_arr; natural_frames_sem];
    
end

% unnatural images
unnatural_frames_mean_arr = [];
unnatural_frames_sem_arr = [];

for f=1:frames
    
    reshaped_angles = reshape(angles(:,:,2,:,f),[],1);
    natural_frames_mean = mean(reshaped_angles,'omitnan');
    natural_frames_std = std(reshaped_angles,'omitnan');
    natural_frames_sem = natural_frames_std / sqrt(length(reshaped_angles));
    
    unnatural_frames_mean_arr = [unnatural_frames_mean_arr; natural_frames_mean];
    unnatural_frames_sem_arr = [unnatural_frames_sem_arr; natural_frames_sem];
end
figure;
errorshade([1:frames],natural_frames_mean_arr,natural_frames_sem_arr,'lineprops',{'Color',[0.2 0.2 0.8]});
errorshade([1:frames],unnatural_frames_mean_arr,unnatural_frames_sem_arr,'lineprops',{'Color',[0.8 0.2 0.2]});
legend('natural','unnatural');
title('data');


% shuffled
natural_frames_mean_arr_shuffle = [];
natural_frames_sem_arr_shuffle = [];
shuffled_frames = randperm(frames-1);

for i=1:length(shuffled_frames)
    f = shuffled_frames(i);
    reshaped_angles = reshape(angles(:,:,1,:,f),[],1);
    
    natural_frames_mean = mean(reshaped_angles,'omitnan');
    natural_frames_std = std(reshaped_angles,'omitnan');
    natural_frames_sem = natural_frames_std / sqrt(length(reshaped_angles));
    
    natural_frames_mean_arr_shuffle = [natural_frames_mean_arr_shuffle; natural_frames_mean];
    natural_frames_sem_arr_shuffle = [natural_frames_sem_arr_shuffle; natural_frames_sem];

end

% unnatural images
unnatural_frames_mean_arr_shuffle = [];
unnatural_frames_sem_arr_shuffle = [];
for i=1:length(shuffled_frames)
    f = shuffled_frames(i);    
    reshaped_angles = reshape(angles(:,:,2,:,f),[],1);
    natural_frames_mean = mean(reshaped_angles,'omitnan');
    natural_frames_std = std(reshaped_angles,'omitnan');
    natural_frames_sem = natural_frames_std / sqrt(length(reshaped_angles));
    
    unnatural_frames_mean_arr_shuffle = [unnatural_frames_mean_arr_shuffle; natural_frames_mean];
    unnatural_frames_sem_arr_shuffle = [unnatural_frames_sem_arr_shuffle; natural_frames_sem];
end
figure;
errorshade([1:frames-1],natural_frames_mean_arr_shuffle,natural_frames_sem_arr_shuffle,'lineprops',{'Color',[0.2 0.2 0.8]});
errorshade([1:frames-1],unnatural_frames_mean_arr_shuffle,unnatural_frames_sem_arr_shuffle,'lineprops',{'Color',[0.8 0.2 0.2]});
legend('natural','unnatural');
title('shuffled');
xlabel('frames');
ylabel('cos(f,pc1)');

