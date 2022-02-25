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

%%


% calculate noise_covariance and mean_responses of each combinations
for initial=1:frames
    
    typeofzooms = 2;
    categories = 3;
    movies = 10;
    frames = 11;
    noise_covariance_pca = num2cell(nan(num_recordings,typeofzooms,categories,movies,frames));
    mean_responses = num2cell(nan(num_recordings,typeofzooms,categories,movies,frames));
    pcs = nan(num_recordings,typeofzooms,categories,movies,frames);
    angles = nan(num_recordings,typeofzooms,categories,movies,frames);
    shuffled_angles = nan(num_recordings,typeofzooms,categories,movies,frames);
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
            cell_res_allframe = removeNan(neuralresponse,z,c,m,frames);
            
            mean_responses = [];
            % get noise covariance matrix, pca, and mean_responses
            for f=1:frames
                cell_res = cell_res_allframe(f,:,:);
                cell_res = squeeze(cell_res);
                [coeff,score,latent,~,explained] = pca(cell_res');
                noise_covariance_pca{i,z,c,m,f} = coeff;
             % mean responses
                f_stimulus = mean(cell_res,2,'omitnan');
                mean_responses = [mean_responses f_stimulus];
            end
            
            if(~isnan(mean_responses))
            for f=1:frames
                shuffled_index= randi([1 frames]);

                noise = noise_covariance_pca{i,z,c,m,f};
                pc = noise(:,1);
                
                noise_ = noise_covariance_pca{i,z,c,m,initial};
                pc_ = noise_(:,1);
                
                
                noise_shuffled = noise_covariance_pca{i,z,c,m,shuffled_index};
                pc_shuffled = noise_shuffled(:,1);
    
                
                angle = pc'*pc_;
                shuffled_angle = pc'*pc_shuffled;
                
                angles(i,z,c,m,f)=angle;
                shuffled_angles(i,z,c,m,f)=shuffled_angle;

            end
            end
    end
    end
    end
end
    drawDataandShuffled(angles,shuffled_angles,frames);

end

