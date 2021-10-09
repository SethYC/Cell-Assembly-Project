% Automated Extraction of cell assembly data Part 1
% 
% Part 1 of the anlysis pipeline for the Russo et al. CA detection results on the
% Lesley recording dataset. Searches for and extracts relevant information
% into a table that is saved and used in the part 2 script. Analysis
% is to answer the following two questions:
%   -task 1 = testing group differences on distribution of cell assemblies
%   (CA) across temporal bins
%   -task 2 = testing group differences on distribution of cell assemblies
%   across number of member neurons per CA
% 
% Currently the script is set up to extract unpruned results from the pyramidal 
% neuron only files, but can be tweaked to look at interneuron/pyramidal
% files, and/or the pruned results by uncommenting and commenting the
% appropiate lines.
%
% Seth Campbell - Aug 1, 2021

% tic %times the script, took roughly a few hours to run through 7191 files

% cd 'Y:\Jennifer\eyeblink_neuron_subset\30_neurons\'  %for pyramidal and interneurons results
cd 'Y:\Jennifer\eyeblink_neuron_subset_pyramidal\30_neurons\' %for pyramidal neurons only results

%file datastore searches the current directory for all the distance type pruned result files
%note: @matfile is used to load the files as i only need two variables from
%each file, thus dramatically saving memory and loading time
ds = fileDatastore("*\*\main_assembly.mat","ReadFcn",@matfile,"FileExtensions",".mat"); 
ds_size = length(ds.Files); 

table_size = ds_size*9; %number of file times 9 temporal bins for each to become number of rows in table

%init a table with preallocated memory, used for storing results 
rats = cell(table_size,1); %e.g. initialize a cell array that has ds_size rows, and 1 column
day = cell(table_size,1);
shuffle = cell(table_size,1);
temporal_bin = cell(table_size,1);

rest1_counts = cell(table_size,1);
task1_counts = cell(table_size,1);
rest2_counts = cell(table_size,1);
rest3_counts = cell(table_size,1);
task2_counts = cell(table_size,1);
rest4_counts = cell(table_size,1);

bin_names = [0.003,0.005,0.01,0.015,0.025,0.04,0.06,0.08,0.10];

t = table(shuffle,rats,day,temporal_bin,rest1_counts,task1_counts,rest2_counts,rest3_counts,task2_counts,rest4_counts); %create a table with these 5 columns
tic
for i = 1:ds_size  %iterate over each file in the datastore    
    i %was for showing progress  
    workspace = read(ds); %read one file from the datastore, when called in the next iteration it gets the next file 
    parameters = matfile([ds.Folders{1}, '\parameters.mat']); %epoch time info is in another file, but always the same folder as main_assembly.mat
    epochs = parameters.epochs; %get only epochs variable out
    epochs_hist_edges = [epochs.Rest1 epochs.Maze1 epochs.Rest2 epochs.Rest3 epochs.Maze2 epochs.Rest4];
    
    %note:because i'm using matfile() to load only the variables i want (instead of the entire workspace), i want to import the Amatrix variable once, 
    %then reuse it, or else if i recall workspace.Amtrix each time i need a value from Amatrix, it slows the script down a ton!
    assembly = workspace.assembly;

    %extract rat id and day from file name
    [filePath] = fileparts(char(ds.Files(i))); %fileparts() extracts the filename as its 2nd output, and char() is applied to convert a cell array to a string
    rat_id = filePath(end-9:end-6);
    day_num = filePath(end-1:end);
    shuffle_num = filePath(end-12:end-11);
    
    %loop through all temporal bins
    for j = 1:length(assembly.bin)
        %init epoch count vars to act as running totals (i.e. r1 = rest 1)
        r1 = 0; t1 = 0; r2 = 0; r3 = 0; t2 = 0; r4 = 0; 
        
        bin = assembly.bin{j}; %current temporal bin
        if isempty(bin) %sometimes there will be no CA found within a bin, so move to next bin
            t((i-1)*9+j,:) = {{shuffle_num}, {rat_id}, {day_num}, {bin_names(j)}, {r1}, {r2}, {t1}, {r3}, {t2}, {r4}}; %set results to zero
            continue
        end
        
        %loop through all CA's in the temporal bin and find activation times
        for k = 1:length(bin.n)
            CA = bin.n{k}; %current CA
        
            %find CA activation timestamps
            %   Note: this is meat of the script: finds the POSITION of whenever a CA 
            %   fires in CA.Time, then uses that to index into bin_edges to get the relative 
            %   timestamp! although a CA could fire at any time WITHIN the
            %   bin specified by two adjacent elements of bin_edges, we
            %   choose to assign the left edge as the time of activation as
            %   a simple approximation.
            activation_ts = bin.bin_edges(find(CA.Time == 1)); 
            activation_ts = activation_ts*1e4; %convert to sec (times 10 to msec, then times 1000 to sec)
            %find counts of CA activations within each of the 6 epochs
            epochs_hist = histcounts(activation_ts, epochs_hist_edges);
            r1 = r1 + epochs_hist(1);
            t1 = t1 + epochs_hist(3);
            r2 = r2 + epochs_hist(5);
            r3 = r3 + epochs_hist(7);
            t2 = t2 + epochs_hist(9);
            r4 = r4 + epochs_hist(11);
        end 
        
        %TO ADD: clolumns for new data (6 of them)
        t((i-1)*9+j,:) = {{shuffle_num}, {rat_id}, {day_num}, {bin_names(j)}, {r1}, {r2}, {t1}, {r3}, {t2}, {r4}}; %init a row of the table with the data we extracted
    end      
end
toc
%convert shuffle number, rats and days columns to doubles instead of a cell array of chars,
%first by iterating through each row of the first three columns of the table
shuffle_double = zeros(1,length(t{:,1}));
rats_double = zeros(1,length(t{:,1}));
day_double = zeros(1,length(t{:,1}));
for i = 1:length(t{:,1})
    foo = t{i,1};
    foo = str2double(char(foo{1}));
    shuffle_double(i) = foo; %create a new array of days as doubles
    
    foo = t{i,2}; %extract a 1x1 cell array for the rat id entry
    foo = str2double(char(foo{1})); %convert to a double, seems convoluted and it is, probably could improve this line or understand better...
    rats_double(i) = foo; %create a new array of rat ids as doubles
    
    foo = t{i,3};
    foo = str2double(char(foo{1}));
    day_double(i) = foo; %create a new array of days as doubles 
end

%create a new table to overwrite the old one
%note: first three columns vars are transposed to become columns as input to table(), and 
%the other arguments are just the third and fourth column of the original table  
% t = table(shuffle_double',rats_double',day_double',t.temporal_bin,t.rest1_counts,t.task1_counts,t.rest2_counts,t.rest3_counts,t.task2_counts,t.rest4_counts); 
t = removevars(t,{'shuffle' 'rats' 'day'});
t = addvars(t,shuffle_double',rats_double',day_double','Before','temporal_bin');

t.Properties.VariableNames = {'shuffle #' 'rat' 'day' 'temporal bin' 'rest 1' 'task 1', 'rest 2', 'rest 3' 'task 2' 'rest 4'}; %name the table columns

%change directories and save the table variable there
cd 'C:\Users\seth.campbell\OneDrive - University of Lethbridge\Documents\Masters\Comp Neuro Project\scripts\'
save('CA counts per epochs unpruned pyramidal only','t')    

% toc