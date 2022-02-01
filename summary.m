close all;
clear;

addpath(genpath('./library'));

%% conditions
NUM_SCALE               = 2;
NUM_CATEGORIES          = 3;
NUM_MOVIES              = 10;
NUM_FRAMES              = 11;
NUM_POPULATIONS         = 9;

STIM_IMAGE_ON_SEC       = 0.2;
STIM_IMAGE_OFF_SEC      = 0.1;

%% enumerations
NATURAL_TYPE_ENUM       = 1;
SYNTHETIC_TYPE_ENUM     = 2;
CONTRAST_TYPE_ENUM      = 3; % N/A


%% LOAD DATA
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

%% CONSTRUCT TABLES
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
for i = 1:num_recordings
    uniq_labels{i}      = char(uniq_label_cat(i));
end

num_table_entries       = numel(recordingLabel);

curvature_model         = nan(num_table_entries, 1);
curvature_model_pixel   = nan(num_table_entries, 1);
relStraighteningModel   = nan(num_table_entries, 1);
distancesModel          = nan(num_table_entries, 1);
corr_Data_vs_V1model    = nan(num_table_entries, 1);

corr_Data2fold          = nan(num_table_entries, 1);

recordingID             = nan(num_table_entries, 1);
num_units               = nan(num_table_entries, 1);
ditherID                = nan(num_table_entries, 1);
scaleID                 = nan(num_table_entries, 1);
categoryID              = nan(num_table_entries, 1);
movieID                 = nan(num_table_entries, 1);
movieLabel              = cell(num_table_entries, 1);

mean_unit_counts        = cell(num_table_entries, 1); % mean spike counts
mean_diff_unit_counts   = cell(num_table_entries, 1); % mean diff b/w frames
median_unit_counts      = cell(num_table_entries, 1);
sigma_g_values          = cell(num_table_entries, 1);
noise_correlations      = nan(num_table_entries, 1);

mean_weights_sd              = nan(num_table_entries, 1); % debugging

for i = 1:num_recordings
    
    % load spike data for each table entry
    i_recording_label           = strtok(recordings(i).name, '.');
    i_recording_label_mask      = i_recording_label == recordingLabel;
    i_filepath                  = fullfile(recordings(i).folder, recordings(i).name);
    load(i_filepath, 'naturalStruct');
    
    for category = categories
        for dither = 1:2
        
            % index from table
            label_mask          = uniq_labels{i} == recordingLabel;
            category_mask       = category == category_list;
            dither_mask         = dither == dither_list;
            i_index             = label_mask & category_mask & dither_mask;
            i_num_movies        = sum(i_index);     % no. movies across zoom1x and zoom2x (e.g. 10 x 2)
            num_orig_movies     = i_num_movies/2;   % no. movies from original source (e.g. 10)
            
            % sort model results to match table labels
            monkey_table        = uniq_labels{i}(1:3);
            penetration_table   = uniq_labels{i}(5:8);
            
            recording_idx       = find(strcmp(uniq_labels{i}(1:3), monkey_id_list) & strcmp(uniq_labels{i}(5:8), penetration_id_list));
            
            recordingID(i_index)= recording_idx * ones(i_num_movies,1);
            num_units(i_index)  = num_units_recording(i);
            
            movie_idx           = (1:num_orig_movies)';
            i_movieID           = [];
            for s = 1:NUM_SCALE
                i_movieID       = cat(1, i_movieID, movie_idx + (s-1)*num_orig_movies);
            end
            movieID(i_index)    = i_movieID;
            if(category == NATURAL_TYPE_ENUM)
                movieLabel(i_index) = natural_movie_labels(i_movieID);
            elseif(category == SYNTHETIC_TYPE_ENUM)
                movieLabel(i_index) = artificial_movie_labels(i_movieID);
            elseif(category == CONTRAST_TYPE_ENUM)
                movieLabel(i_index) = contrast_movie_labels(i_movieID);
            end
            
            % associate spike data to each table entry
            i_recording_index           = find(i_index);
            
            mean_weights_sd(i_recording_index)   = repmat(weights_sd(i), [numel(i_recording_index), 1]);
            
            % mean (and median) spike counts for each frame
            movie_x_units_mean_counts    = [];
            movie_x_units_median_counts  = [];
            for s = 1:NUM_SCALE
                for m = 1:num_orig_movies
                    index_location      = m + (s-1)*num_orig_movies;
                    
                    means                       = squeeze(nanmean(  naturalStruct.sortedResponses(s, category, m, :, :, :), repeats_dim));
                    mean_unit_counts{i_recording_index(index_location)}     = means;
                    
                    medians                     = squeeze(nanmedian(naturalStruct.sortedResponses(s, category, m, :, :, :), repeats_dim));
                    median_unit_counts{i_recording_index(index_location)}   = medians;
                    
                end
            end
            
            % mean delta spike-counts (b/w frames)
            mean_diff       = nan(i_num_movies, naturalStruct.num_units);
            for s = 1:NUM_SCALE
                for m = 1:num_orig_movies
                    tmp_mean_resp           = [];
                    for d = 1:dither
                        movie_mean_resp     = squeeze(nanmean(naturalStruct.sortedResponses(s, category, m, d:dither:end, :, :), repeats_dim));
                        tmp_mean_resp       = cat(1, tmp_mean_resp, diff(movie_mean_resp,1));
                       
                    end
                    
                    index_location      = m + (s-1)*num_orig_movies;
                    mean_diff_unit_counts{i_recording_index(index_location)}     = tmp_mean_resp;
                    
                end
            end
            
            % noise correlations
            movie_ns    = nan(i_num_movies, 1);
            for s = 1:NUM_SCALE
                for m = 1:num_orig_movies
                    d_ns    = nan(dither, 1);
                    for d = 1:dither
                        movie_resp      = squeeze(naturalStruct.sortedResponses(s, category, m, d:dither:end, :, :));
                        movie_ns_matrix = getNoiseCorrelation_singleMovie(movie_resp);
                        upper_tri       = triu(movie_ns_matrix, 1); 
                        d_ns(d)         = nanmean(upper_tri(:));
                    end
                    movie_ns(m + (s-1)*num_orig_movies)     = nanmean(d_ns);
                end
            end
            noise_correlations(i_recording_index)   = movie_ns;
            
            
            % book keeping
            ditherID(i_index)   = dither * ones(i_num_movies, 1);
            categoryID(i_index) = category * ones(i_num_movies, 1);
            
            i_scaleID           = [];
            for s = 1:NUM_SCALE
                i_scaleID       = cat(1, i_scaleID, s * ones(num_orig_movies, 1));
            end
            scaleID(i_index)    = i_scaleID;
            
            i_curv_model        = [];
            for s = 1:NUM_SCALE
                i_curv_model    = cat(1, i_curv_model, squeeze(c_V1_model(recording_idx,dither,s,category,:)));
            end
            curvature_model(i_index)    = i_curv_model;
            
            i_curv_pixel        = [];
            for s = 1:NUM_SCALE
                i_curv_pixel    = cat(1, i_curv_pixel, squeeze(c_V1_pixel(recording_idx,dither,s,category,:)));
            end
            curvature_model_pixel(i_index)    = i_curv_pixel;
            
            i_delta_curv_model  = [];
            for s = 1:NUM_SCALE
                i_delta_curv_model  = cat(1, i_delta_curv_model, squeeze(delta_curvature(recording_idx,dither,s,category,:)));
            end
            relStraighteningModel(i_index)    = i_delta_curv_model;
            
            i_distances         = [];
            for s = 1:NUM_SCALE
                i_distances     = cat(1, i_distances, squeeze(distance_model(recording_idx,dither,s,category,:)));
            end
            distancesModel(i_index)    = i_distances;
            
            corr_model          = [];
            for s = 1:NUM_SCALE
                corr_model      = cat(1, corr_model, squeeze(corr_V1_model_data(recording_idx,dither,s,category,:)));
            end
            corr_Data_vs_V1model(i_index) = corr_model;
            
            corr_data           = [];
            for s = 1:NUM_SCALE
                corr_data       = cat(1, corr_data, squeeze(data_correlations.corr_pearson(s,category,:,recording_idx)));
            end
            corr_Data2fold(i_index)    = corr_data;
            
            
        end
    end
end

% generate table for each sequence type
% label/re-order things to enable logical masking by each variable name
% (summary_table.Properties.VariableNames')
summary_table           = table;

% data ID (TODO: indicate movie ID)
summary_table           = addvars(summary_table, recordingLabel);
summary_table           = addvars(summary_table, recordingID);
summary_table           = addvars(summary_table, num_units);
summary_table           = addvars(summary_table, ditherID);
summary_table           = addvars(summary_table, scaleID);
summary_table           = addvars(summary_table, categoryID);
summary_table           = addvars(summary_table, movieID);
summary_table           = addvars(summary_table, movieLabel);

% neural data
deltaCurvatureNeural    = Acutedatacombined235x.DeltaCurvatureneural;
deltaCurvatureControl   = Acutedatacombined235x.DeltaCurvaturecontrol;
relStraighteningNeural  = deltaCurvatureNeural - deltaCurvatureControl;
distancesNeural         = Acutedatacombined235x.Distances;
corr_Data_vs_Pred       = Acutedatacombined235x.Ratescorrelationempiricalvspredicted;
corr_Data_vs_Pred_var   = Acutedatacombined235x.Variancecorrelationempiricalvspredicted;
corr_Data_vs_Pred_covar = Acutedatacombined235x.Covariancecorrelationempiricalvspredicted;
summary_table           = addvars(summary_table, deltaCurvatureNeural);
summary_table           = addvars(summary_table, deltaCurvatureControl);
summary_table           = addvars(summary_table, relStraighteningNeural);
summary_table           = addvars(summary_table, distancesNeural);
pixelCurvatures         = Acutedatacombined235x.Pixelcurvature;
summary_table           = addvars(summary_table, pixelCurvatures);
summary_table           = addvars(summary_table, corr_Data_vs_Pred);
summary_table           = addvars(summary_table, corr_Data_vs_Pred_var);
summary_table           = addvars(summary_table, corr_Data_vs_Pred_covar);
summary_table           = addvars(summary_table, mean_unit_counts);
summary_table           = addvars(summary_table, mean_diff_unit_counts);
summary_table           = addvars(summary_table, median_unit_counts);
summary_table           = addvars(summary_table, noise_correlations);

% data quality
summary_table           = addvars(summary_table, corr_Data2fold);

% image computable model
summary_table           = addvars(summary_table, curvature_model);
summary_table           = addvars(summary_table, curvature_model_pixel);
summary_table           = addvars(summary_table, relStraighteningModel);
summary_table           = addvars(summary_table, distancesModel);
summary_table           = addvars(summary_table, corr_Data_vs_V1model);

summary_table           = addvars(summary_table, mean_weights_sd);

bootstraps_neural       = bootstrap_table.bootstraps_neural;
bootstraps_control      = bootstrap_table.bootstraps_control;

summary_table           = addvars(summary_table, bootstraps_neural);
summary_table           = addvars(summary_table, bootstraps_control);

save('summary_table.mat', 'summary_table', ...
    'num_units_recording', ...
    'natural_movie_labels', 'artificial_movie_labels');
