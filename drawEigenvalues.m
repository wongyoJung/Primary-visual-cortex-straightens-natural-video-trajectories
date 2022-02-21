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
%% draw eigenvalues

typeofzooms = 2;
categories = 3;
movies = 10;
frames = 11;
% calculate noise_covariance and mean_responses of each combinations
num_of_units = nan(num_recordings,typeofzooms,categories,movies);
for i = 1:num_recordings
    
    % load spike data for each table entry
    i_recording_label           = strtok(recordings(i).name, '.');
    i_recording_label_mask      = i_recording_label == recordingLabel;
    i_filepath                  = fullfile(recordings(i).folder, recordings(i).name);
    load(i_filepath, 'naturalStruct');
    
    neuralresponse = naturalStruct.sortedResponses;    
    eigval_natural = [];
    eigval_unnatural = [];
    min_units = 0;

    for z=1:typeofzooms  
    for c=1:2
    for m=1:movies

            nanrows = [];
            nancolumns = [];
            % remove repeats having nan
            for f=1:frames
                cell_resps = squeeze(neuralresponse(z,c,m,f,:,:));
                [nanrow,nancolumn] = find(isnan(cell_resps));
                   
                nanrow = unique(nanrow);  
                nancolumn = unique(nancolumn);  

                nanrows = [nanrows; nanrow];       
                nancolumns = [nancolumns; nancolumn];              

            end
            nanrows= unique(nanrows);
            cell_res_allframe = squeeze(neuralresponse(z,c,m,:,:,:));
            cell_res_allframe(:,nanrows,:)=[];
            units = size(cell_res_allframe,2);
            
            if(min_units==0)
                min_units = units;
            else
                if(min_units>units)
                    min_units = units;
                end
            end
            
            eigvals= [];
            % get noise covariance matrix, pca, and mean_responses
            for f=1:frames
                cell_res = cell_res_allframe(f,:,:);
                cell_res = squeeze(cell_res);
                [coeff,score,latent,~,explained] = pca(cell_res');
                eigvals = [eigvals latent];
                eigvals = eigvals./sum(eigvals);
            end  
            
            if(c==1)
                eigval_natural = cat(3, eigval_natural, eigvals(1:min_units,:));
            else
                eigval_unnatural = cat(3, eigval_unnatural,  eigvals(1:min_units,:));
            end

    end
    end
    end

figure; 

mean_val = squeeze(mean(eigval_natural,[2 3]));
std_val = squeeze(std(eigval_natural,[],[2 3]));
subplot(2,1,1);
errorshade([1:min_units], mean_val,std_val, 'lineprops',{'Color',[0.2 0.2 0.8]});
set(gca,'yscale','log');

title('natural');
xlabel("principal component");
ylabel("Eigenvalues");
xlim([1 min_units]);


mean_val = squeeze(mean(eigval_unnatural,[2 3]));
std_val = squeeze(std(eigval_unnatural,[],[2 3]));

subplot(2,1,2);
errorshade([1:min_units], mean_val,std_val,'lineprops',{'Color',[0.8 0.2 0.2]});
set(gca,'yscale','log');
title('unnatural');
xlabel("principal component");
ylabel("Eigenvalues");
xlim([1 min_units]);



end


