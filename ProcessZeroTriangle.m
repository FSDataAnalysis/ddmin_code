t%% vector of dmin distance ddm
%% vector of distance amplitude dAW distance_dWA
%% vector of height Amplitude height_dWA

% clear all;
% load mydata;


prompt= {'cutt off ddmin:'};
dlg_title='Inputs';
num_lines=1;
default={'0.5e-9'};
answer=inputdlg(prompt,dlg_title,num_lines,default);

%%%%%%% Data on outliers/or start and end of vectors %%%%%%%%%%%%%%%%%%%%%%


cut_off_ddmin=str2double(answer(1));  % this is the 

foo_ddmin=ddm;

remove_points=(foo_ddmin(1,:))>=cut_off_ddmin;

counter_removals=sum(remove_points);

new_ddm=ddm(remove_points);
new_distance_dAW=distance_dWA(remove_points);
new_height_dAW=height_dWA(remove_points);
