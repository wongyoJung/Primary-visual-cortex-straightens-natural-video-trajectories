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


%% pca

typeofzooms = 2;
categories = 3;
movies = 10;
frames = 11;
noise_covariance_pca = num2cell(nan(num_recordings,typeofzooms,categories,movies,frames));
mean_responses = num2cell(nan(num_recordings,typeofzooms,categories,movies,frames));
angles = nan(num_recordings,typeofzooms,categories,movies,frames);
% calculate noise_covariance and mean_responses of each combinations
num_of_units = nan(num_recordings,typeofzooms,categories,movies);
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
%             cell_res_allframe(:,:,nancolumns)=[];


            % get noise covariance matrix, pca, and mean_responses
            for f=1:frames
                cell_res = cell_res_allframe(f,:,:);
                cell_res = squeeze(cell_res);
                [coeff,score,latent,~,explained] = pca(cell_res');
                noise_covariance_pca{i,z,c,m,f} = coeff;
             % mean responses
                f_stimulus = mean(cell_res,2,'omitnan');
                mean_responses{i,z,c,m,f} = f_stimulus;
            end
            
            
            if(~isnan(mean_responses{i,z,c,m,f}))
            for f=1:frames-1
                f_diff = mean_responses{i,z,c,m,f+1} - mean_responses{i,z,c,m,f};
                norm_f_diff = f_diff./norm(f_diff);
                
                noise = noise_covariance_pca{i,z,c,m,f};
                pc = noise(:,3);
%                 proj = pc'*norm_f_diff;
                angle = pc'*norm_f_diff;
%                 angle = norm(proj)^2;
                angles(i,z,c,m,f)=angle;
            end
            end
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
    x = zeros(size(reshaped_angles))+f;
%     hold on; scatter(x,reshaped_angles);  
%     hold on; scatter([f],[mean(num_of_units)]);
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
xlabel('frames');
ylabel('cos(f,pc3)');

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




%% replace nan with average of each unit's mean(spike count)
typeofzooms = 2;
categories = 3;
movies = 10;
frames = 11;
noise_covariance_pca = num2cell(nan(num_recordings,typeofzooms,categories,movies,frames));
mean_responses = num2cell(nan(num_recordings,typeofzooms,categories,movies,frames));
angles = nan(num_recordings,typeofzooms,categories,movies,frames);
% calculate noise_covariance and mean_responses of each combinations
num_of_units = nan(num_recordings,typeofzooms,categories,movies);
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
            % get noise covariance matrix, pca, and mean_responses
            for f=1:frames
                cell_res = cell_res_allframe(f,:,:);
                cell_res = squeeze(cell_res);
                [coeff,score,latent,~,explained] = pca(cell_res');
                noise_covariance_pca{i,z,c,m,f} = coeff;
             % mean responses
                f_stimulus = mean(cell_res,2,'omitnan');
                mean_responses{i,z,c,m,f} = f_stimulus;
            end
            
            
            if(~isnan(mean_responses{i,z,c,m,f}))
            for f=1:frames-1
                f_diff = mean_responses{i,z,c,m,f+1} - mean_responses{i,z,c,m,f};
                norm_f_diff = f_diff./norm(f_diff);
                
                noise = noise_covariance_pca{i,z,c,m,f};
                pc = noise(:,1);
%                 proj = pc'*norm_f_diff;
                angle = pc'*norm_f_diff;
%                 angle = norm(proj)^2;
                angles(i,z,c,m,f)=angle;
            end
            end
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

%% alignment between f_diff, based on geometric location, and pc of noise covariance

typeofzooms = 2;
categories = 3;
movies = 10;
frames = 11;
noise_covariance_pca = num2cell(nan(num_recordings,typeofzooms,categories,movies,frames));
mean_responses = num2cell(nan(num_recordings,typeofzooms,categories,movies,frames));
pcs = nan(num_recordings,typeofzooms,categories,movies,frames);
angles = nan(num_recordings,typeofzooms,categories,movies,frames);
shuffled_angles = nan(num_recordings,typeofzooms,categories,movies,frames);

% calculate noise_covariance and mean_responses of each combinations

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
                f_mean = mean_responses(:,f);
                f_mean_ = repmat(f_mean,[1, frames]);
                f_diff = mean_responses - f_mean_;
                norm_f_diff = vecnorm(f_diff);
                [out,idx]=sort(norm_f_diff);
                min_index=idx(2);
                
                f_adjacent = mean_responses(:,f)-mean_responses(:,min_index);
                noise = noise_covariance_pca{i,z,c,m,f};
                pc = noise(:,1);
                
                list=setdiff(1:frames,f);
                shuffled_index=list(randperm(frames-1,1));
                f_shuffle = mean_responses(:,f)-mean_responses(:,shuffled_index);
                
                
                angle = pc'*f_adjacent;
                shuffled_angle = pc'*f_shuffle;
                
                angles(i,z,c,m,f)=angle;
                shuffled_angles(i,z,c,m,f)=shuffled_angle;
                
                
                
            end
            end
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

for f=1:frames

    shuffle_reshaped_angles = reshape(shuffled_angles(:,:,1,:,f),[],1);
    natural_frames_mean = mean(shuffle_reshaped_angles,'omitnan')
    natural_frames_std = std(shuffle_reshaped_angles,'omitnan');
    natural_frames_sem = natural_frames_std / sqrt(length(shuffle_reshaped_angles));
    
    natural_frames_mean_arr_shuffle = [natural_frames_mean_arr_shuffle; natural_frames_mean];
    natural_frames_sem_arr_shuffle = [natural_frames_sem_arr_shuffle; natural_frames_sem];

end

% unnatural images
unnatural_frames_mean_arr_shuffle = [];
unnatural_frames_sem_arr_shuffle = [];
for f=1:frames
    reshaped_angles = reshape(shuffled_angles(:,:,2,:,f),[],1);
    natural_frames_mean = mean(reshaped_angles,'omitnan');
    natural_frames_std = std(reshaped_angles,'omitnan');
    natural_frames_sem = natural_frames_std / sqrt(length(reshaped_angles));
    
    unnatural_frames_mean_arr_shuffle = [unnatural_frames_mean_arr_shuffle; natural_frames_mean];
    unnatural_frames_sem_arr_shuffle = [unnatural_frames_sem_arr_shuffle; natural_frames_sem];
end
figure;
errorshade([1:frames],natural_frames_mean_arr_shuffle,natural_frames_sem_arr_shuffle,'lineprops',{'Color',[0.2 0.2 0.8]});
errorshade([1:frames],unnatural_frames_mean_arr_shuffle,unnatural_frames_sem_arr_shuffle,'lineprops',{'Color',[0.8 0.2 0.2]});
legend('natural','unnatural');
title('shuffled');
xlabel('frames');
ylabel('cos(f,pc1)');