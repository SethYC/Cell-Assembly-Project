% load('CA counts per epochs unpruned pyramidal only (Aug 2)');
load('CA probability per epochs unpruned pyramidal only (Aug 2)');


bin_names = [0.003,0.005,0.01,0.015,0.025,0.04,0.06,0.08,0.10];

%convert temporal bin and later columns to numbers and not cell arrays
a = cell2mat(t.("temporal bin"));
b = cell2mat(t.("CA #")); %uncomment later
c = cell2mat(t.("rest 1"));
d = cell2mat(t.("task 1"));
e = cell2mat(t.("rest 2"));
f = cell2mat(t.("rest 3"));
g = cell2mat(t.("task 2"));
h = cell2mat(t.("rest 4"));

t = removevars(t,{'temporal bin' 'CA #' 'rest 1' 'task 1', 'rest 2', 'rest 3' 'task 2' 'rest 4'}); %uncomment later
% t = removevars(t,{'temporal bin' 'rest 1' 'task 1', 'rest 2', 'rest 3' 'task 2' 'rest 4'});

t = addvars(t,a,b,c,d,e,f,g,h,'After','day'); %uncomment later
% t = addvars(t,a,c,d,e,f,g,h,'After','day');

t.Properties.VariableNames = {'shuffle #' 'rat' 'day' 'temporal bin' 'CA #' 'rest 1' 'task 1', 'rest 2', 'rest 3' 'task 2' 'rest 4'}; %uncomment later
% t.Properties.VariableNames = {'shuffle #' 'rat' 'day' 'temporal bin' 'rest 1' 'task 1', 'rest 2', 'rest 3' 'task 2' 'rest 4'};


%extracts only a certain temporal bin, comment out if you want full results
% idx = t.("temporal bin") == bin_names(9);
% t(~idx,:) = [];

%collapse rows across 9 bins for each shuffle/rat/day combo
[groups, tid] = findgroups(t(:,1:3));
r1 = splitapply(@sum,t.("rest 1"),groups);
t1 = splitapply(@sum,t.("task 1"),groups);
r2 = splitapply(@sum,t.("rest 2"),groups);
r3 = splitapply(@sum,t.("rest 3"),groups);
t2 = splitapply(@sum,t.("task 2"),groups);
r4 = splitapply(@sum,t.("rest 4"),groups);

t = addvars(tid, r1, t1, r2, r3, t2, r4); %create new table from results
t.Properties.VariableNames = {'shuffle #' 'rat' 'day' 'rest 1' 'task 1', 'rest 2', 'rest 3' 'task 2' 'rest 4'};

% t2 = []; %init empty array/table
% for i = 1:9:length(t.rat)
%     t_temp = t(i:i+8,:); 
%     t2 = [t2;t_temp];
% end
    
%reshape matrix so that shuffles is a 3rd dimensions


%find mean across third dimension (across shuffles)


%~~~
[groups, tid] = findgroups(t(:,2:3));
r1 = splitapply(@mean,t.("rest 1"),groups);
t1 = splitapply(@mean,t.("task 1"),groups);
r2 = splitapply(@mean,t.("rest 2"),groups);
r3 = splitapply(@mean,t.("rest 3"),groups);
t2 = splitapply(@mean,t.("task 2"),groups);
r4 = splitapply(@mean,t.("rest 4"),groups);

%include the following if you want a binary version of the table
% r1(r1>=1) = 1;
% t1(t1>=1) = 1;
% r2(r2>=1) = 1;
% r3(r3>=1) = 1;
% t2(t2>=1) = 1;
% r4(r4>=1) = 1;

%add custom vars
all_rest = (r1+r2+r3+r4)/4;
all_task = (t1+t2)/2;
pre_task_rest = (r1+r3)/2;
post_task_rest = (r2+r4)/2;

t = addvars(tid, r1, t1, r2, r3, t2, r4, all_rest, all_task, pre_task_rest, post_task_rest); %create new table from results
t.Properties.VariableNames = {'rat' 'day' 'rest 1' 'task 1', 'rest 2', 'rest 3' 'task 2' 'rest 4' 'all rest' 'all task' 'pre-task rest' 'post-task rest'};

% writetable(t,'collapsed CA probability per epochs unpruned pyramidal only (Aug 2).csv');

%create 6 epoch tables

%init 11x34 matrix (rats as rows, days as columns + 1 column of rat numbers)
r1_chart = NaN(11,34);
t1_chart = NaN(11,34);
r2_chart = NaN(11,34);
r3_chart = NaN(11,34);
t2_chart = NaN(11,34);
r4_chart = NaN(11,34);
ar_chart = NaN(11,34);
at_chart = NaN(11,34);
pretr_chart = NaN(11,34);
posttr_chart = NaN(11,34);

rat_ids = [ 8417,
            8419,
            8564,
            8645,
            8646,
            8700,
            8778,
            8820,
            8886,
            8957,
            8981,];

for i = 1:length(t.rat)
    id = t.rat(i);
    day = t.day(i);
    
    rat_pos = find(rat_ids == id);
    r1_chart(rat_pos,day+1) = t.("rest 1")(i);
    t1_chart(rat_pos,day+1) = t.("task 1")(i);
    r2_chart(rat_pos,day+1) = t.("rest 2")(i);
    r3_chart(rat_pos,day+1) = t.("rest 3")(i);
    t2_chart(rat_pos,day+1) = t.("task 2")(i);
    r4_chart(rat_pos,day+1) = t.("rest 4")(i);
    ar_chart(rat_pos,day+1) = t.("all rest")(i);
    at_chart(rat_pos,day+1) = t.("all task")(i);
    pretr_chart(rat_pos,day+1) = t.("pre-task rest")(i);
    posttr_chart(rat_pos,day+1) = t.("post-task rest")(i);
end

r1_chart(:,1) = rat_ids;
t1_chart(:,1) = rat_ids;
r2_chart(:,1) = rat_ids;
r3_chart(:,1) = rat_ids;
t2_chart(:,1) = rat_ids;
r4_chart(:,1) = rat_ids;
ar_chart(:,1) = rat_ids;
at_chart(:,1) = rat_ids;
pretr_chart(:,1) = rat_ids;
posttr_chart(:,1) = rat_ids;

col_names = {'rat #' '1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' '21' '22' '23' '24' '25' '26' '27' '28' '29' '30' '31' '32' '33'};

r1_table = array2table(r1_chart,'VariableNames',col_names);
t1_table = array2table(t1_chart,'VariableNames',col_names);
r2_table = array2table(r2_chart,'VariableNames',col_names);
r3_table = array2table(r3_chart,'VariableNames',col_names);
t2_table = array2table(t2_chart,'VariableNames',col_names);
r4_table = array2table(r4_chart,'VariableNames',col_names);
ar_table = array2table(ar_chart,'VariableNames',col_names);
at_table = array2table(at_chart,'VariableNames',col_names);
pretr_table = array2table(pretr_chart,'VariableNames',col_names);
posttr_table = array2table(posttr_chart,'VariableNames',col_names);

writetable(r1_table,'CA probability per epochs unpruned pyramidal - rest 1 table (Aug 2).csv');
writetable(t1_table,'CA probability per epochs unpruned pyramidal - task 1 table (Aug 2).csv');
writetable(r2_table,'CA probability per epochs unpruned pyramidal - rest 2 table (Aug 2).csv');
writetable(r3_table,'CA probability per epochs unpruned pyramidal - rest 3 table (Aug 2).csv');
writetable(t2_table,'CA probability per epochs unpruned pyramidal - task 2 table (Aug 2).csv');
writetable(r4_table,'CA probability per epochs unpruned pyramidal - rest 4 table (Aug 2).csv');
writetable(ar_table,'CA probability per epochs unpruned pyramidal - all rest table (Aug 2).csv');
writetable(at_table,'CA probability per epochs unpruned pyramidal - all task table (Aug 2).csv');
writetable(pretr_table,'CA probability per epochs unpruned pyramidal - pre-task rest table (Aug 2).csv');
writetable(posttr_table,'CA probability per epochs unpruned pyramidal - post-task rest table (Aug 2).csv');

% writetable(r1_table,'CA count per epochs unpruned pyramidal - rest 1 table (Aug 2).csv');
% writetable(t1_table,'CA count per epochs unpruned pyramidal - task 1 table (Aug 2).csv');
% writetable(r2_table,'CA count per epochs unpruned pyramidal - rest 2 table (Aug 2).csv');
% writetable(r3_table,'CA count per epochs unpruned pyramidal - rest 3 table (Aug 2).csv');
% writetable(t2_table,'CA count per epochs unpruned pyramidal - task 2 table (Aug 2).csv');
% writetable(r4_table,'CA count per epochs unpruned pyramidal - rest 4 table (Aug 2).csv');
% writetable(ar_table,'CA count per epochs unpruned pyramidal - all rest table (Aug 2).csv');
% writetable(at_table,'CA count per epochs unpruned pyramidal - all task table (Aug 2).csv');
% writetable(pretr_table,'CA count per epochs unpruned pyramidal - pre-task rest table (Aug 2).csv');
% writetable(posttr_table,'CA count per epochs unpruned pyramidal - post-task rest table (Aug 2).csv');

%how to specifiy a certain temporal bin for the table
% idx = t.("temporal bin") == 0.0100
% t2 = t(idx,:);