% Automated Extraction of cell assembly data Part 2
% 
% Part 2 of the anlysis pipeline for the Russo et al. CA detection results on the
% Lesley recording dataset. Using the table of data from part 1, gathers the information
% creates charts, and runs stats to help answer the two questions below. 
%   -task 1 = testing group differences on distribution of cell assemblies
%   (CA) across temporal bins
%   -task 2 = testing group differences on distribution of cell assemblies
%   across number of member neurons per CA
%
% Seth Campbell - July 10, 2021

%import table from part 1 of anlaysis (see Automated_Cell_Assembly_Results_Extraction_Part_1.m)
load('30 neuron controlled pyramidal only output table (July 6 all).mat')

%define rat groups based on their id #, e.g. young vs old, high NCA vs. low NCA, etc.
young = [8419,8570,8645,8646,8820,8957];
old = [8417,8564,8700,8778,8886,8981];
high_NCA = [8417,8419,8564,8645,8646,8957]; 
low_NCA = [8570,8700,8778,8820,8886,8981];
remap = [8417, 8419, 8700, 8820, 8886, 8981];
non_remap = [8564, 8570, 8645, 8646, 8778, 8957];

bins = [0.003,0.005,0.01,0.015,0.025,0.04,0.06,0.08,0.10]; %the 9 preset bins of temporal precision, for task 1
edges = 0.5:9.5; %edges for the histogram  x-axes for task 2's charts

%% compute results based on groupings, specifically young vs. old

[task_1_results_young_unsummed, task_2_results_young] = groupTally(t,young); %groupTally is defined at the bottom (with other local functions)
[task_1_results_old_unsummed, task_2_results_old] = groupTally(t,old);

%create a histogram of the results for task 2
ax1 = subplot(3,2,1); %note: using subplot here so i can have multiple figures in one window to compare, also note that all task 2 charts will be on the left side in the final plot
histogram(task_2_results_young, edges, 'Normalization', 'probability') %histogram from the collected member neuron counts, but tallied such that its normalized (*see note at end on normalization)
hold on;
histogram(task_2_results_old, edges, "Normalization", 'probability');
title("Histogram of member neuron counts per cell assembly")
xlabel("member neurons per cell assembly")
ylabel("Perentage of total counts")
legend("Young", "Old")
hold off

%create a histogram of the results for task 1
ax2 = subplot(3,2,2); %note: task 1 charts will always be in the right side in the end plot
a = histcounts(task_1_results_young_unsummed,bins); %unlike task 2, a bar chart needs to be used for task 1 due to the discrete x-axis, this normalizes the results manually
bar(a/sum(a),"FaceAlpha",.6); %gonna overlay two bar charts, so changing the alpha helps them visually overlap
hold on;
a = histcounts(task_1_results_old_unsummed,bins);
bar(a/sum(a),"FaceAlpha",.6);
title("Histogram of Cell Assemblies per Bin Size");
xlabel("Bin Sizes");
set(gca,'xticklabel',{0.003,0.005,0.01,0.015,0.025,0.04,0.06,0.08,0.10}) %set label for each bar
ylabel("Perentage of total counts");
legend('Young', 'Old','Location','northwest')
hold off;

%% high vs low NCA results

[task_1_results_high_NCA_unsummed, task_2_results_high_NCA] = groupTally(t,high_NCA);
[task_1_results_low_NCA_unsummed, task_2_results_low_NCA] = groupTally(t,low_NCA);

%create and label a histogram of the results for task 2
ax3 = subplot(3,2,3);
histogram(task_2_results_high_NCA, 'Normalization', 'probability')
hold on
histogram(task_2_results_low_NCA, 'Normalization', 'probability')
xlabel("member neurons per cell assembly")
ylabel("Perentage of total counts")
legend("High Cell Assemblies", "Low Cell Assemblies")
hold off

%create and label a histogram of the results for task 1
ax4 = subplot(3,2,4);
a = histcounts(task_1_results_high_NCA_unsummed,bins);
bar(a/sum(a),"FaceAlpha",.6);
hold on;
a = histcounts(task_1_results_low_NCA_unsummed,bins);
bar(a/sum(a),"FaceAlpha",.6);
xlabel("Bin Sizes");
set(gca,'xticklabel',{0.003,0.005,0.01,0.015,0.025,0.04,0.06,0.08,0.10})
ylabel("Perentage of total counts");
legend('High Cell Assemblies', 'Low Cell Assemblies', 'Location','northwest')
hold off;

%% remapping vs non-remapping results

[task_1_results_non_remap_unsummed, task_2_results_non_remap] = groupTally(t,non_remap);
[task_1_results_remap_unsummed, task_2_results_remap] = groupTally(t,remap);

%create and label a histogram of the results for task 2
ax5 = subplot(3,2,5);
histogram(task_2_results_non_remap, 'Normalization', 'probability')
hold on
histogram(task_2_results_remap, 'Normalization', 'probability')
xlabel("member neurons per cell assembly")
ylabel("Perentage of total counts")
legend("Non-remapping", "Remapping")
hold off

%create and label a histogram of the results for task 1
ax6 = subplot(3,2,6);
a = histcounts(task_1_results_non_remap_unsummed,bins);
bar(a/sum(a),"FaceAlpha",.6);
hold on;
a = histcounts(task_1_results_remap_unsummed,bins);
bar(a/sum(a),"FaceAlpha",.6);
xlabel("Bin Sizes");
set(gca,'xticklabel',{0.003,0.005,0.01,0.015,0.025,0.04,0.06,0.08,0.10})
ylabel("Perentage of total counts");
legend('Non-remapping', 'Remapping', 'Location','northwest')
hold off;

linkaxes([ax1,ax3,ax5]); %match the y-axes of the charts on the left
linkaxes([ax2,ax4,ax6]); %match the y-axes of the charts on the right
sgtitle("30 Neuron Controlled Unpruned Cell Assemblies") %title of the entire subplot

%% stats

[p,h,stats] = ranksum(task_2_results_young,task_2_results_old); %Wilcoxon Rank Sum test between two groups
fprintf('Young vs. Old member neurons: p = %f, zval = %f\n', p, stats.zval); %display results in command window

[p,h,stats] = ranksum(task_2_results_high_NCA,task_2_results_low_NCA);
fprintf('High Cell Assemblies vs. Low Cell Assemblies member neurons: p = %f, zval = %f\n', p, stats.zval);

[p,h,stats] = ranksum(task_2_results_non_remap,task_2_results_remap);
fprintf('Non-remapping vs. Remapping member neurons: p = %f, zval = %f\n', p, stats.zval);

[p,h,stats] = ranksum(task_1_results_young_unsummed,task_1_results_old_unsummed);
fprintf('Young vs. Old temporal bins: p = %f, zval = %f\n', p, stats.zval);

[p,h,stats] = ranksum(task_1_results_high_NCA_unsummed,task_1_results_low_NCA_unsummed);
fprintf('High Cell Assemblies vs. Low Cell Assemblies temporal bins: p = %f, zval = %f\n', p, stats.zval);

[p,h,stats] = ranksum(task_1_results_non_remap_unsummed,task_1_results_remap_unsummed);
fprintf('Non-remapping vs. Remapping temporal bins: p = %f, zval = %f\n', p, stats.zval);

%% functions

% Iterates over a given group of rat's data and collects the info from the table 
% to answer task 1 and task 2's questions 
% 
% note: commented out code is deprecated, and was used with tally1(), the 
% function used to output task_1_results instead of task_1_results_unsummed
%
% input: 
%   -table "t" of collected rat data from a previous script (e.g.
%   June_task_version.m)
%   -vector "group" containing the rat ids for a group such as young or old,
%   etc. 
% output: 
%   -task_1_results_unsummed: a long vector where each element is the # neurons found in a cell
%   assembly, for all cell assemblies for the group
%   -task_2_results: also a long vector where each element is the # neurons found in a cell
%   assembly, for all cell assemblies for the group
function [task_1_results_unsummed, task_2_results] = groupTally(t,group)
%     task_1_results = [0 0 0 0 0 0 0 0 0];
    task_1_results_unsummed = [];
    task_2_results = [];
%     task_1_total = length(group); %init number of rats in group, and will be adjusted if any rat is found missing
%     j=1;
    for i = group %iterate over each rat in the group
%         result = tally1(t,i);
%         if isempty(result)
%             task_1_total = task_1_total - 1; %if a rat was not present in the results table, then when averaging later, the total # of rats to divide by needs to be adjusted!
%         else
%             task_1_results = task_1_results + result;
%         end
        task_1_results_unsummed = [task_1_results_unsummed, tally1_unsummed(t,i)]; %build up a vector for each rat in the group
        task_2_results = [task_2_results, tally2(t,i)];
%         j = j + 1;
    end      

%     task_1_results = task_1_results/task_1_total; %averaging results by number of rats found
end

%Aids task 2, returns the total count of neurons per CA for a single rat over all its days
%
% input: 
%   -table "t" of collected rat data
%   -rat # "rat_id" to specifiy which rat
% output:
%   -row vector of all CA member count entries for the given rat
function output_task_2 = tally2(t,rat_id)
    idx = t.rat == rat_id; %create a logical vector of row positions in the table where rat #'s matches rat_id
    t2 = t(idx,:);  %use the logical vector to filter/extract those relevant results into a new table
    results = t2{:,4}; %extract results from the 4th column (neurons_per_CA)
    output_task_2 = [results{:}]; %collapse results into 1 dimension
end

%**(DEPRECATED)** Aids task 1 by tallying the occurence of cell assemblies in each bin size
%for a given rat (determined by "rat_id", extracted from table "t")
% function result = tally1(t,rat_id)
%     bins = [0.003,0.005,0.01,0.015,0.025,0.04,0.06,0.08,0.10];
%     
%     idx = t.rat == rat_id;
%     t2 = t(idx,:); 
%     
%     if isempty(t2.rat) % handle cases when a certain rat is not found (e.g. none of its recordings met a new criteria)
%         result = [];
%         return
%     end
%     
%     results = t2.CA_per_bin;
%     output_task_1 = [results{:}];
%     output_task_1_categorical = categorical(output_task_1,bins);
%     result = histcounts(output_task_1_categorical);
% end

%Supercedes the role of tally1(), operating the same way as tally2(), but
%with temporal bin information instead of member neuron counts
%
% input: 
%   -table "t" of collected rat data
%   -rat # "rat_id" to specifiy which rat
% output:
%   -row vector of all CA_per_bin entries for the given rat
function output_task_1_unsummed = tally1_unsummed(t,rat_id)
    idx = t.rat == rat_id;
    t2 = t(idx,:); 
    results = t2{:,5}; %extract results from the 5th column (CA_per_bin)
    output_task_1_unsummed = [results{:}];
end

%%
% *Normalization Note*:
%   -original results did not use normalization, but to control for the
%   influence of recording files with more CA on the results, we used a
%   probabiliy normalization method, such that the height of all bars for a
%   group in a chart sums to 1, thus measuring relative precence in each
%   bin/column. 
%   -another version of this entire script exists that first applied this
%   normalization procedure on the individual recording files, THEN
%   averaged the results. This version produced weaker group differences 
%   across all charts, but is technically a better anlaysis in that it 
%   controlled for the influence of number of CA found per recording file 
%   before averaging results across all recording files (because we found
%   that recording files with more total CA's found also had more large
%   CA's for example). It was not used for the SFN poster however because
%   we were unsure of how to properly apply stats on this version. 