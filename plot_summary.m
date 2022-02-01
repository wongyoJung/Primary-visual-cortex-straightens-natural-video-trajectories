close all;
clear;

addpath(genpath('./library'));

%% load summarized data from 'summary.m'
load('summary_table.mat', 'summary_table','num_units_recording',...
        'natural_movie_labels', 'artificial_movie_labels', 'contrast_movie_labels');
%%%%%%%%%%%%%% table variables %%%%%%%%%%%%%%
%%%% ID's %%%%
% summary_table.recordingLabel
% summary_table.recordingID
% summary_table.ditherID
% summary_table.categoryID
% summary_table.scaleID
% summary_table.movieID
%
%%%% neural variables %%%%
% summary_table.deltaCurvatureNeural
% summary_table.deltaCurvatureControl
% summary_table.relStraighteningNeural
% summary_table.distancesNeural
%
%%%% data quality %%%%
% summary_table.corrData2fold
%
%%%% image computable model %%%%
% summary_table.relStraighteningModel
% summary_table.distancesModel
% summary_table.corrData_vs_V1model

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
%% PLOTTING CONVENTIONS
hsv_colors              = hsv;
recording_colors        = hsv_colors(1 : floor(size(hsv, 1) / numel(9)) : size(hsv, 1), :);
movie_colors            = hsv_colors(1 : floor(size(hsv, 1) / NUM_MOVIES) : size(hsv, 1), :);

red                     = [255 0 0 ]  ./256;
orange                  = [255 128 0] ./ 256;
green                   = [0 165 80]  ./ 256;
blue                    = [0 175 235] ./ 256;
purple                  = [145 40 145]./ 256;   
grey                    = [150 150 150] ./ 256;

myColorMap              = zeros(201, 3);
myColorMap(1, :)        = [.5 .5 .5];
myColorMap(2:end, 1)    = [linspace(0,1,100) ,      ones(1, 100)];
myColorMap(2:end, 2)    = [linspace(0,1,100) , linspace(1,0,100)];
myColorMap(2:end, 3)    = [     ones(1, 100) , linspace(1,0,100)];

fontname                = 'Helvetica';
fontsize                = 20;
fontangle               = 'oblique';

scale_labels            = {'zoom1x', 'zoom2x'};
cat_labels              = {'natural', 'artificial', 'contrast'};


%% feature variables
categories              = [NATURAL_TYPE_ENUM, SYNTHETIC_TYPE_ENUM];
category_labels         = {'NATURAL', 'ARTIFICIAL'};
category_colors         = {blue, green};
% main results from dither = 2 
dither_main              = 2;

category_list           = summary_table.categoryID;
dither_list             = summary_table.ditherID;
scale_list              = summary_table.scaleID;
population_list         = summary_table.recordingID;
movie_list              = summary_table.movieID;

bootstraps_neural       = summary_table.bootstraps_neural;
bootstraps_control      = summary_table.bootstraps_control;

scale_id                = unique(scale_list)';
scale_id(isnan(scale_id))=[];
num_scales              = numel(scale_id);

movie_id                = unique(summary_table.movieID);
invalid_movies          = isnan(movie_id);
movie_id(invalid_movies)=[];

num_movies              = numel(movie_id);
num_zoom1x_movies       = num_movies/2;

population_id           = unique(population_list);
% remove recordings w/ hardware failure
population_id(2)        = [];
population_id(7)        = [];
population_id(isnan(population_id))=[];
num_populations         = numel(population_id);

%% apply data inclusion criteria
DISTANCE_THRESH         = 0.25;
MEAN_PREDICTION_THRESH  = 0.75;
VAR_PREDICTION_THRESH   = 0.5;
COVAR_PREDICTION_THRESH = 0.25;

distance_filter         = summary_table.distancesNeural      > DISTANCE_THRESH;
var_expl_filter_mean    = summary_table.corr_Data_vs_Pred.^2 > MEAN_PREDICTION_THRESH;
var_expl_filter_var     = summary_table.corr_Data_vs_Pred_var.^2 > VAR_PREDICTION_THRESH;
var_expl_filter_covar   = summary_table.corr_Data_vs_Pred_covar.^2 > COVAR_PREDICTION_THRESH;

filtered_index_mean     = distance_filter & var_expl_filter_mean;
filtered_index_mean_var = distance_filter & var_expl_filter_mean & var_expl_filter_var;
filtered_index          = distance_filter & var_expl_filter_mean & var_expl_filter_var & var_expl_filter_covar;

% fix exclusion to naturals & dither = 2
dither2                 = (category_list == NATURAL_TYPE_ENUM) & (dither_list == 2);
dither1                 = (category_list == NATURAL_TYPE_ENUM) & (dither_list == 1);


%% stimulus
load('stimulus/stim_info.mat', ...
    'natural_movie_labels', 'artificial_movie_labels', 'contrast_movie_labels', ...
    'natural_movie_frame', 'artificial_movie_frame', 'contrast_movie_frame');


%% Summary matrices
summary_matrices        = nan(numel(categories), numel(dither_main), num_movies , num_populations);
distance_matrices       = nan(numel(categories), numel(dither_main), num_movies , num_scales, num_populations);
dCurv_neural            = nan(numel(categories), numel(dither_main), num_movies , num_populations);
dCurv_control           = nan(numel(categories), numel(dither_main), num_movies , num_populations);
model_matrices          = nan(numel(categories), numel(dither_main), num_movies , num_populations);

isNeuralOutsideNull95   = nan(numel(categories), numel(dither_main), num_movies , num_populations);

NATURAL_MAX_DEG         = 75;
ARTIFICIAL_MAX_DEG      = 75;

histogram_max_deg       = 180;
histogram_num_bins      = 36;

figure_folder           = fullfile('./figures/');


for d = dither_main
    
    d_mask          = dither_list == d;
    
    for c = categories
        
        fig_summary = figure('units', 'normalized', 'outerposition', [0,0,0.9,0.9]);
    
        if(c == 1)
            max_degrees = NATURAL_MAX_DEG;
        else
            max_degrees = ARTIFICIAL_MAX_DEG;
        end
        
        c_mask      = category_list == c;
        
        index_mask      = c_mask & d_mask & filtered_index;
        
        summary_matrix  = nan(num_movies , num_populations);
        
        % zoom1x 
        zoom1x_mask                     = index_mask & (scale_list == scale_id(1));
        for m = 1:num_movies
            movie_zoom1x_mask           = zoom1x_mask & (movie_list == movie_id(m));
            for p = 1:num_populations
                pop_mask                = movie_zoom1x_mask & (population_list == population_id(p));
                if(sum(pop_mask)>0)
                    summary_matrix(m,p)     = summary_table.relStraighteningNeural(pop_mask);
                    dCurv_neural(c,d,m,p)   = summary_table.deltaCurvatureNeural(pop_mask);
                    dCurv_control(c,d,m,p)  = summary_table.deltaCurvatureControl(pop_mask);
                    model_matrices(c,d,m,p) = summary_table.relStraighteningModel(pop_mask);
                    
                    distance_matrices(c,d,m,1,p) = summary_table.distancesNeural(pop_mask);
                    
                    mu_neural               = mean(bootstraps_neural(pop_mask,:));
                    range95                 = quantile(bootstraps_control(pop_mask,:), [0.0275, 0.975]);
                    if( c == 1)
                        isTrue  = (mu_neural < range95(1));
                    else
                        isTrue  = (mu_neural > range95(end));
                    end
                    isNeuralOutsideNull95(c,d,m,p)  = isTrue;
                    
                end
            end
        end
        
        % zoom2x
        zoom2x_mask                     = index_mask & (scale_list == scale_id(2));
        for m = 1:num_movies
            movie_zoom2x_mask           = zoom2x_mask & (movie_list == movie_id(m));
            for p = 1:num_populations
                pop_mask                = movie_zoom2x_mask & (population_list == population_id(p));
                if(sum(pop_mask)>0)
                    summary_matrix(m,p)     = summary_table.relStraighteningNeural(pop_mask);
                    dCurv_neural(c,d,m,p)   = summary_table.deltaCurvatureNeural(pop_mask);
                    dCurv_control(c,d,m,p)  = summary_table.deltaCurvatureControl(pop_mask);
                    model_matrices(c,d,m,p) = summary_table.relStraighteningModel(pop_mask);

                    distance_matrices(c,d,m,1,p) = summary_table.distancesNeural(pop_mask);
                    
                    mu_neural               = mean(bootstraps_neural(pop_mask,:));
                    range95                 = quantile(bootstraps_control(pop_mask,:), [0.0275, 0.975]);
                    if( c == 1)
                        isTrue  = (mu_neural < range95(1));
                    else
                        isTrue  = (mu_neural > range95(end));
                    end
                    isNeuralOutsideNull95(c,d,m,p)  = isTrue;
                    
                end
            end
        end
        
        
        % sort by natural movies only
        if( c == 1 )
            [sorted_pop, sorted_pop_idx]    = sort(nanmean(summary_matrix,1), 'descend');
            [sorted_mov, sorted_mov_idx]    = sort(nanmean(summary_matrix,2), 'descend');
        elseif(c == 2) 
            marg_pop    = nanmean(summary_matrix,1);
            marg_mov    = nanmean(summary_matrix,2);
            
            sorted_pop  = marg_pop(sorted_pop_idx);
            sorted_mov  = marg_mov(sorted_mov_idx);
            
            % keep track of most/least entangling in artificial sequences
            [arti_sorted_pop, arti_sorted_pop_idx]    = sort(nanmean(summary_matrix,1), 'ascend');
            [arti_sorted_mov, arti_sorted_mov_idx]    = sort(nanmean(summary_matrix,2), 'ascend');
        end
        
        summary_matrices(c,d,:,:)       = summary_matrix;
        
        sorted_summary_matrix           = summary_matrix(:,sorted_pop_idx);
        sorted_summary_matrix           = sorted_summary_matrix(sorted_mov_idx,:);
        new_pop_idx                     = 1:numel(sorted_pop_idx);
        new_mov_idx                     = 1:numel(sorted_mov_idx);
        
        valid_entries                   = ~isnan(sorted_summary_matrix);
        se_pop                          = sqrt(nanvar(sorted_summary_matrix, [], 1) ./ sum(valid_entries, 1));
        se_mov                          = sqrt(nanvar(sorted_summary_matrix, [], 2) ./ sum(valid_entries, 2));

        % summary matrix
        axes('position', [-0.05 + 0.2*(c-1), 0.1, 0.35, 0.85]);
        title(sprintf('%s (step-size: %d) \n %d/%d entries (%2.0f%% yield)', ...
            category_labels{c}, d, ...
            sum(~isnan(summary_matrix(:))), numel(summary_matrix(:)), 100*sum(~isnan(summary_matrix(:)))/numel(summary_matrix(:))));
        xlabel(sprintf('population ID'));
        ylabel(sprintf('movie ID'));
        phyplot([], [], 'k-', ...
            'yticks',       1:size(sorted_summary_matrix,1), ...
            'yticklabels',  sprintfc('%d', flipud(new_mov_idx)), ...
            'xticks',       1:size(sorted_summary_matrix,2), ...
            'xticklabels', sprintfc('%d', new_pop_idx), ...
            'linewidth', 1, 'width', 1, 'fontsize', 12);
        hold on;
        
        
        %max_degrees = 10+max(round(abs(sorted_summary_matrix(:))));
        clipped_sorted_summary_matrix = sorted_summary_matrix;
        clip_idx    = abs(sorted_summary_matrix) > max_degrees;
        clipped_sorted_summary_matrix(clip_idx) = (max_degrees-1) .* sign(sorted_summary_matrix(clip_idx));
        imagesc(flipud(clipped_sorted_summary_matrix)); 
        axis equal;
        colormap(myColorMap);
        caxis([-max_degrees max_degrees]);

        marker_size = 11;
        
        % population marginals
        axes('position', [0.45, 0.1, 0.3, 0.3]);
        
        if(c == 1)
            title(sprintf('Populations (marginals)'));
            ylabel(sprintf('straightening (deg)'));
            xlabel(sprintf('population ID'));
        end
        phyplot([], [], 'k-', ...
            'yticks',       [-max_degrees,0,max_degrees], ...
            'xticks',       new_pop_idx, ...
            'xticklabels',  sprintfc('%d', new_pop_idx), ...
            'linewidth', 1, 'width', 1, 'fontsize', 12);
        hold on;
        for i = 1:numel(new_pop_idx)
            text(i, max_degrees*1.1, sprintf('n=%d',num_units_recording(sorted_pop_idx(i))));
        end
        
        plot(1:numel(sorted_pop), 0 * ones(1,numel(sorted_pop)), '--', 'color', 'b');
        verrorbar(1:numel(sorted_pop), sorted_pop, se_pop, '-', 'color', 'k', 'linewidth', 1.5);
        plot(new_pop_idx, sorted_pop, 'wo','markerfacecolor', 'k', 'markersize', marker_size, 'linewidth', 1.5);
        
       
        % movie marginals
        axes('position', [0.25, 0.5, 0.7, 0.4]);
        if(c == 1)
            title(sprintf('Movies (marginals)'));
            ylabel(sprintf('straightening (deg)'));
            xlabel(sprintf('movie ID'));
        end
        phyplot([], [], 'k-', ...
            'yticks',       [-max_degrees,0,max_degrees], ...
            'xticks',       new_mov_idx, ...
            'xticklabels',  sprintfc('%d', new_mov_idx), ...
            'linewidth', 1, 'width', 1, 'fontsize', 12);
        hold on;
        
        plot(1:numel(sorted_mov), 0 * ones(1,numel(sorted_mov)), '--', 'color', 'b');
        verrorbar(1:numel(sorted_mov), sorted_mov, se_mov, '-', 'color', 0.5*ones(1,3), 'linewidth', 1.5);
        plot(new_mov_idx, sorted_mov, 'o','markerfacecolor', 'w', 'markersize', marker_size,  'linewidth', 1.5, 'color', 0.5*ones(1,3));
        
        if(c == 1)
            plot(new_mov_idx(1),              sorted_mov(1),              'wo','markerfacecolor', 0.5*ones(1,3), 'markersize', marker_size, 'linewidth', 1.5);
            plot(new_mov_idx(num_movies/2+1), sorted_mov(num_movies/2+1), 'wo','markerfacecolor', 'r', 'markersize', marker_size, 'linewidth', 1.5);
            plot(new_mov_idx(end),            sorted_mov(end),            'wo','markerfacecolor', 'k', 'markersize', marker_size, 'linewidth', 1.5);
        else
            plot(arti_sorted_mov_idx(1),              arti_sorted_mov(1),              'wo','markerfacecolor', 0.5*ones(1,3), 'markersize', marker_size, 'linewidth', 1.5);
            plot(arti_sorted_mov_idx(num_movies/2+1), arti_sorted_mov(num_movies/2), 'wo','markerfacecolor', 'r', 'markersize', marker_size, 'linewidth', 1.5);
            plot(arti_sorted_mov_idx(end),            arti_sorted_mov(end),            'wo','markerfacecolor', 'k', 'markersize', marker_size, 'linewidth', 1.5);
        end
        
        % print movie labels for most/medium/least straightening movie
        if(c == NATURAL_TYPE_ENUM)
            sorted_movie_labels = natural_movie_labels(sorted_mov_idx);
            fprintf('\n------DITHER %d, %s -------------\n', d, category_labels{c});
            fprintf('\nLEAST   straightening movie: %s\n', sorted_movie_labels{1});
            fprintf('\nMEDIUM straightening movie: %s\n', sorted_movie_labels{num_movies/2+1});
            fprintf('\nMOST  straightening movie: %s\n', sorted_movie_labels{num_movies});
        elseif(c == SYNTHETIC_TYPE_ENUM)
            sorted_movie_labels = artificial_movie_labels(arti_sorted_mov_idx);
            fprintf('\n------DITHER %d, %s -------------\n', d, category_labels{c});
            fprintf('\nLEAST   entangling movie: %s\n', sorted_movie_labels{1});
            fprintf('\nMEDIUM entangling movie: %s\n', sorted_movie_labels{num_movies/2});
            fprintf('\nMOST  entangling movie: %s\n', sorted_movie_labels{num_movies});
        end
              
        % print out whether empirical estimate is plausible under null
        % distribution
        tmp             = isNeuralOutsideNull95(c,d,:,:);
        tmp             = tmp(:);
        tmp(isnan(tmp)) = [];
        perc_outside_95 = 100 * sum(tmp)/numel(tmp);
        fprintf('\n%s-DITHER %d\n%2.1f %% of sequences fell outside the 95%% interval of null-distribution\n(%d out of %d)\n', ...
            category_labels{c}, d, perc_outside_95, sum(tmp), numel(tmp));
        
        
        % histograms
        axes('position', [0.75, 0.55 - 0.5*(c-1),  0.25, 0.35]);
        num_bins            = histogram_num_bins;
        rel_straight        = summary_matrix(:);
        rel_straight(isnan(rel_straight)) = [];
        [p_wilcoxon, h]     = signrank(rel_straight);
        
        p_limit     = 0.001;
        if(p_wilcoxon > p_limit)
            title(sprintf('%s (median: %2.2f^o)\n p_{wilcoxon}: %1.3f, N=%d', category_labels{c}, median(rel_straight), p_wilcoxon, numel(rel_straight)));
        else
            title(sprintf('%s (median: %2.2f^o)\n p_{wilcoxon} < %1.3f, N=%d', category_labels{c}, median(rel_straight), p_limit, numel(rel_straight)));
        end
        xlabel(sprintf('\\Delta curvature (deg)'));
        ylabel(sprintf('Counts'));
        
        % bins
        max_deg                 = histogram_max_deg;
        edges                   = linspace(-max_deg, max_deg, num_bins);
        tmp                     = diff(edges);
        bin_width               = tmp(1);
        bin_centers             = linspace(edges(1), edges(end), num_bins-1);
        
        % counts
        max_y                   = 50;
        yticks                  = [0  max_y];
        yticklabels             = {num2str(0), num2str(max_y)};
        
        phyplot([], [], 'k-', ...
            'xticks',       [-max_deg, -max_deg/2, 0, max_deg/2, max_deg ], ...
            'yticks',       yticks, ...
            'yticklabels',  yticklabels, ...
            'linewidth', 1, 'width', 1, 'fontsize', 14);
        
        hold on;
        plot([0,0], [0, max_y], 'k--', 'linewidth', 1);
        counts_V1           = histcounts(rel_straight, edges);
        b                   = bar(bin_centers, counts_V1, 'FaceColor', category_colors{c}, 'EdgeColor', 'k', 'LineWidth', 1.5);
        plot(nanmedian(rel_straight),  max_y, 'v', 'markerfacecolor', category_colors{c}, 'markeredgecolor', 'k', 'markersize', 20);
          
        filepath    = fullfile(figure_folder, sprintf('summary_%s.pdf', cat_labels{c}));
        orient landscape;
        set(fig_summary, 'PaperSize', [22 14]);
        print(fig_summary, filepath, '-dpdf', '-r0');
    
    end
    
end

