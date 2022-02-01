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


%% =============================== noise correlation within a single unit ===============================

%     'Response tensor for image sequences:                  '
%     'dimension 1: 2 scales (zoom1x, zoom2x)                '
%     'dimension 2: 3 category (natural, synthetic, contrast)'
%     'dimension 3: 10 movies                                '
%     'dimension 4: 11 frames per movie                      '
%     'dimension 5: sorted units                             '
%     'dimension 6: repetitions                              '

typeofzooms = 2;
categories = 3;
movies = 10;
frames = 11;
noise_correlation = nan(num_recordings,typeofzooms,categories,movies,frames);



for i = 1:num_recordings
    
    % load spike data for each table entry
    i_recording_label           = strtok(recordings(i).name, '.');
    i_recording_label_mask      = i_recording_label == recordingLabel;
    i_filepath                  = fullfile(recordings(i).folder, recordings(i).name);
    load(i_filepath, 'naturalStruct');
    
    neuralresponse = naturalStruct.sortedResponses;    

    for z=1:typeofzooms
    for c=1:categories
    for m=1:movies
            for f=1:frames
                        cell_resp =  neuralresponse(z,c,m,f,:,:);
                        cell_resps = squeeze(cell_resp);
                        cond_mean       = mean(cell_resps,2,'omitnan');
                        cond_sd         = mean(cell_resps,2,'omitnan');
                        residual        = cell_resps - cond_mean;
                        z_score        = residual./cond_sd;
                        mean_zscores = mean(z_score,'all','omitnan');
                        noise_correlation(i,z,c,m,f) = mean_zscores;
            end
    end
    end
    end
end

%% draw graph
figure;
for c=1:categories
    categorized_resp = reshape(noise_correlation(:,:,c,:,:),[],11);
    meann = mean(categorized_resp,1,'omitnan');
    stdd = std(categorized_resp,1,'omitnan');
    semm = stdd./sqrt(size(categorized_resp,1));
    hold on;
    if(c==1)
    errorshade([1:11],meann, semm,'lineProps',{'Color',[0.2 0.2 0.8]});
    elseif(c==2)
    errorshade([1:11],meann, semm,'lineProps',{'Color',[0.8, 0.2, 0.2]});

    end
end
    legend('natural','unnatural','contrast');
    xlabel("frames");
    ylabel("z-score");
    title(" noise correlation within a single unit ")

%% =============================== noise correlation between pairs of neurons ===============================


cross_noise_correlation = nan(num_recordings,typeofzooms,categories,movies,frames);

for i = 1:num_recordings
    
    % load spike data for each table entry
    i_recording_label           = strtok(recordings(i).name, '.');
    i_recording_label_mask      = i_recording_label == recordingLabel;
    i_filepath                  = fullfile(recordings(i).folder, recordings(i).name);
    load(i_filepath, 'naturalStruct');
    
    neuralresponse = naturalStruct.sortedResponses;    

    for z=1:typeofzooms
    for c=1:categories
    for m=1:movies
            for f=1:frames
                      cell_resps = squeeze(neuralresponse(z,c,m,f,:,:));
                      corrcoef_pair = corrcoef(cell_resps');
                      colvect= corrcoef_pair(find(triu(corrcoef_pair,1)));
                      corr   = mean(colvect,'omitnan');
                      cross_noise_correlation(i,z,c,m,f)=corr;
            end
    end
    end
    end
end


figure;
% natural movies
drawNoiseCorrChangeByTime(cross_noise_correlation,1,frames,[0.2 0.2 0.8])
% unnatural movies
drawNoiseCorrChangeByTime(cross_noise_correlation,2,frames,[0.8 0.2 0.2])
legend('natural','unnatural');
xlabel("frames");
ylabel("coefficient");
title(" pearson coefficient between all pair of neurons ")

