%%%%% Add amplitude phase deflection files and it will calculate minimum
%%%%% distance of approach separation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% ddm  This is the vector before removing outliers


%%%% ddm_clean  This is the vector after removing outliers 
%%%% distance_dWA
%%%% height_dWA

clear all
close all
clc
tic

originaldir=pwd;
format shortEng
%The following code should be used for an automated process
% Select the folder containing all of the .txt files in the order you would
% like them to be analyzed- no subfolders!!
fPath = uigetdir('.', 'Select directory containing the files');
if fPath==0, error('no folder selected'), end
fNames = dir(fullfile(fPath,'*.txt') );
fNames = strcat(fPath, filesep, {fNames.name});
count_trials=0;

prompt= {'Smoothing Coefficient:', 'InVolts:', 'Cut off in nm:', 'Remove points start:', 'Remove points end:', ...
    'Remove Outliers (0/1):', 'Cut_off dAW:'};
dlg_title='Inputs';
num_lines=1;
default={'0.03','40', '0.2e-9', '10', '2', '1', '0.5'};
answer=inputdlg(prompt,dlg_title,num_lines,default);

%%%%%%% Data on outliers/or start and end of vectors %%%%%%%%%%%%%%%%%%%%%%

s_d_min=str2double(answer(1)); 
cut_off=str2double(answer(3));  % this is the 
remove_start= str2double(answer(4));
remove_end= str2double(answer(5));
Remove_outliers= str2double(answer(6));
cut_off_dAW= str2double(answer(7));

sub_num=0;

%%%% Smoothing 

s_AmEx=0.04;
s_defl=0.02;
s_d_min_Incr=0.02;


%%%%%%%%%%%%%%%%%%%%%%%% Selection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Smooth= s_d_min;	
AmpInvOLS =  str2double(answer(2));	

M_size=5;						% Marker size. 6 for FFT simulation 3 for continuous
set(0,'DefaultAxesFontSize',14)  %%% This sets the font size in figures
set(0,'DefaultAxesFontSize',14)  %% axis size figures
set(0,'DefaultAxesLinewidth',2)  %% Lne width figures
box on

count_trials=count_trials+1;
%mkdir(num2str(count_trials));

Extension=1;                            % If 1 then it does extension otherwise with 0  retraction curve (experimental only) 
Smooth_raw_data=1;						%%% Smoothens the raw data. Requied 1 for experimental only  This slows down. Also requires to rearrange dmin from larger to smaller. 
Averaging_percentage_of_signal=0.02;    %% This is the percentage of avergaing for setting phase at 90, deflection at zero, etc. DEFAULT 0.1 
crop_initial_data=0;					%%% Allows you to crop both ends of initial data 
Limit_AmplitudeToA0=1;					%%% This limits amplitudes to A0, so there are no complex numbers for zc>>A0 due to thermal noise
Maximum_lenght_data_set=0;				%% If 1 it will chop data to the number especified below 
Maximum_lenght_data=5000;				% If a number is chosen the length of the vectors wil be between that and double that if 1 is chosen above 													
plot_data_processing_figures=1;			%%% If zero you do not see how the data is being prepared to be processed. Figs. 1-2 4-10 not shown and extra analysis 
plot_raw_curves=0;						% If 1 it plots raw curves, deflection, amplitude and phase 
save_files_processing=1;				%%% This saves plotted figures 1-3 in the begining				
	


% Cycle through each file individually


fit_struct=struct;
fitted_line_y=struct;


counter_file=1;
cut_off_is=zeros(1,length(fNames));
D_max_cut_off=zeros(1,length(fNames));
D_min_cut_off=zeros(1,length(fNames));
ddm=zeros(1,length(fNames));

MultiplierError=1;   % this is to find minima in difference of dmin
distance_dWA=[];
height_dWA=[];

for i=1:length(fNames)
    DtAq=load(fNames{i});
%This part of the code makes a new folder for every text file. In this
%version of the code, every time the program saves, it changes its
%directory to that of the new folder, and then changes it back after the
%save
    [pathstr, name, ext]=fileparts(fNames{i});
    foldname=name;
    cd(fPath);
    mkdir(foldname);
    cd(originaldir);
  

% ====================================================================================================================================================
% ===================== USER settings =========================================================================================================
% ====================================================================================================================================================


%%%% More details  ==========================================================================================%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Sorting_zc_larger_to_smaller=0;         % DEFAULT 0. If 1, then it forces zc to always decrease, i.e. drift removal, if 0, it does not do it, i.e. drift present, complex numbers might result. 
Removing_repeat_zc=0;                   % Default 0. IMPORTANT. Sorting Zc required for this to work. Otherwise it removes many points. if 1 then it removes repeated zc.

%%%%%==================================================================================================================================== %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%% Experimental loading detail
%%%%%%%%%%%%%%%% ====================================================================================%%%%%%%%%%%%%%%%%%%%%%%%%%%									  
%%%% standard %%%%
ZEx =1 ; DfEx =2 ; ZRet = 3;  DfRet = 4; AEx =5 ;  ARet =6 ; PEx =7 ;  PRet =8 ;
%%%%%%%%%%%%%%%% END loading details ====================================================================================%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%========= Smoothening coefficients    ==============================================================================================%%%%
MultiplierError=1;  


% ====================================================================================================================================================
% ===================== end of user settings =========================================================================================================
% ====================================================================================================================================================


%%

count_figures=0;
rrr=8;    % Standard 

length_vectors=size(DtAq,2);

for ccc=1:rrr:length_vectors
	
	% Experimental data load------------------------------------------------------------------------------

	set_d=DtAq(:,ccc:ccc+rrr-1);
	
	Pre_processing_raw_data;
	
	
    ZsnrEx =ZsnrEx-ZsnrEx(end);   %% last zc is zero 
	
	
	%%%%%%%%%%%%%% PLOT amplitude, deflection and phase %%%%%%%%%%%%%%%%%%%%%%%%%%
	%close (count_figures+1)
	if plot_data_processing_figures==1
		figure (count_figures+1)
		plot(ZsnrEx,AmEx, '.k', 'Markersize',M_size, 'displayname','Amplitud A');    %%% 'fontsize',14,  ,  , 

		hold on
		title(' Amplitude and deflection versus zc','fontsize',12)
		xlabel('Zc','fontsize',14) 
		ylabel('Amplitude','fontsize',14)
	
		plot(ZsnrEx,DflEx, '.c', 'Markersize',M_size,'displayname','Amplitud A');
        
        if save_files_processing==1
            cd(strcat(fPath,filesep,foldname));
            saveas(count_figures+1, strcat(num2str(count_trials),'_',num2str(1)),'fig')
            cd(originaldir);
        end
	
	end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculate parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
	
	%%%%%%%%%%%%% Conservative branch %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	A0 = mean(AmEx(1:floor(0.02*length(AmEx))));
	AmEx_dummy=AmEx;
	ANGLE_CONSERVATIVE=zeros(length(PhEx),1);
	trial=zeros(length(PhEx),1);
 	for iii=1:length(PhEx)
		
		if AmEx_dummy(iii, 1)>A0    %% fluctuations in A0 are due to thermal noise
			AmEx_dummy(iii, 1)=A0;
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LIMITING MAXIMUM AMPLITUDES TO A0 GETS RID OF COMPLEX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			if Limit_AmplitudeToA0==1
				AmEx(iii,1)=A0;
			end	
		end
	end	
	
	for iii=1:length(PhEx)
		trial(iii,1)=asin(AmEx_dummy(iii, 1)/A0);
		ANGLE_CONSERVATIVE(iii, 1)=trial(iii, 1)*180/pi;
		if PhEx(iii,1)>90
			ANGLE_CONSERVATIVE(iii, 1)=180-ANGLE_CONSERVATIVE(iii, 1);
		end
	end


  %%%%%%%%%%% Phase in zc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if plot_data_processing_figures==1
		figure (count_figures+2)
		plot(ZsnrEx,PhEx, '.r', 'Markersize',M_size, 'displayname','Phase degrees');
		box on
		hold on
		plot(ZsnrEx,ANGLE_CONSERVATIVE, '.b', 'Markersize',M_size, 'displayname','Phase Conservative');
		title(' Phase versus zc','fontsize',12)
		xlabel('Zc','fontsize',14) 
		ylabel('Phase','fontsize',14)
        
        if save_files_processing==1
            cd(strcat(fPath,filesep,foldname));
            saveas(count_figures+2, strcat(num2str(count_trials),'_',num2str(2)),'fig')
            cd(originaldir);
        end
	end
	
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculate parameters II %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		
			
    d_min=ZsnrEx-AmEx;     %%%%%%%  d_min defined raw ====================================
    d_min_Smth= smooth(d_min,s_d_min,'rloess');
    d_min_Incr=d_min(2:end)-d_min(1:end-1);
    d_min_Incr_Smth= smooth(d_min_Incr,s_d_min_Incr,'rloess');

    %%% plot dmin  ====================================================================================%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(count_figures+3)
    hold on
    plot(ZsnrEx, d_min, '.r','Markersize',M_size, 'displayname','Dmin');
    hold on
    box on
    title(' Dmin versus zc','fontsize',12)
    xlabel('Zc','fontsize',14) 
    ylabel('Dmin','fontsize',14)
    plot(ZsnrEx, d_min_Smth, '.b','Markersize',M_size);
    

    %%% END plot dmin  ====================================================================================%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  
    %%% Finding and plotting minimum point in separation  ============================================================================= %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      find_minimum_delta=min(d_min_Smth);
%     element_find_minimum_delta=find((d_min_Smth)==find_minimum_delta);
%     plot(ZsnrEx(element_find_minimum_delta),  d_min_Smth(element_find_minimum_delta),'Vk', 'Markersize',8)  %%% Plots minimum value found  in dmin
%     text(0.7*max(ZsnrEx),0.8*min(d_min_Smth), ['Min Dmin:' num2str(d_min_Smth(element_find_minimum_delta))],'fontsize',12)


    dmin_zero_interaction=d_min_Smth(1:300);
    ZsnrEx_zero_interaction=ZsnrEx(1:300);
    
    
    plot(ZsnrEx_zero_interaction, dmin_zero_interaction, '.r','Markersize',M_size*2);
    
    file_name_foo = sprintf( 'file%i',  counter_file );
     
    [dumb_fit,dumb_gof] = fit(ZsnrEx_zero_interaction,dmin_zero_interaction, 'poly1');
    fit_struct.(file_name_foo)=dumb_fit;
    fitted_line_y.(file_name_foo)=fit_struct.(file_name_foo).p1*ZsnrEx+fit_struct.(file_name_foo).p2;
    
    
    %%%% Linear fit %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot(ZsnrEx, fitted_line_y.(file_name_foo), '-k', 'Markersize',5)
    %%%%% Minimum accross %%%%%
    
    %%%%%  difference from fit to real data %%%%%%%%%
    
    difference_is.(file_name_foo)=abs(d_min_Smth-fitted_line_y.(file_name_foo));
    
    cut_off_is_dumb=find(difference_is.(file_name_foo)>cut_off);

    cut_off_is(counter_file)=cut_off_is_dumb(1);
    
    %%%% Find reduced amplitude Vector %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     position_to_cut=cut_off_is(counter_file); 
%     reduced_ZsnrEx=ZsnrEx(position_to_cut:end);
%     reduced_AmEx=AmEx(position_to_cut:end);
%     reduced_AmEx_Smth= smooth(reduced_AmEx,s_AmEx,'rloess');

    %%% reduced dmin find %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    AmEx_Smth= smooth(AmEx,s_AmEx,'rloess');
    
    reduced_1=cut_off_is(counter_file);
    reduced_AmEx=AmEx_Smth(reduced_1:end);
    reduced_zc=ZsnrEx(reduced_1:end);
    
    reduced_2=floor(cut_off_dAW*length(reduced_AmEx));
    
    reduced_AmEx_2=reduced_AmEx(reduced_2:end);
    reduced_zc_2=reduced_zc(reduced_2:end);
    
    max_is=max(reduced_AmEx_2);
    
    index_max=find(reduced_AmEx_2==max_is);
    
    index_max_is=index_max(1);
    
    distance_dWA(counter_file)=abs(reduced_zc_2(index_max_is)-reduced_zc_2(end));
    height_dWA(counter_file)=reduced_AmEx_2(index_max_is);
    
    figure(count_figures+4)
    hold on
    plot(ZsnrEx, AmEx_Smth, '.k', 'Markersize',5);
    plot(reduced_zc_2(index_max_is), reduced_AmEx_2(index_max_is) , 'Vb', 'Markersize',15);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure(count_figures+3)
    hold on
    D_max_cut_off(counter_file)=d_min_Smth( cut_off_is(counter_file));
    D_min_cut_off(counter_file)=min(d_min_Smth);
    
    ddm(counter_file)=abs(D_max_cut_off(counter_file)-D_min_cut_off(counter_file));
    
    plot(ZsnrEx(cut_off_is(counter_file)-sub_num), d_min_Smth(cut_off_is(counter_file)-sub_num), 'Vk', 'Markersize',10)
    
    
    cd(strcat(fPath,filesep,foldname));
    saveas(count_figures+3, strcat(num2str(count_trials),'_',num2str(3)),'fig');
    saveas(count_figures+4, strcat(num2str(count_trials),'_',num2str(4)),'fig');
    cd(originaldir);
    
end

close all                     % closes all figures

counter_file=counter_file+1;
end


if Remove_outliers==1
    
    Data_input=ddm';

    Chauvenete;
    Outliers_ddm_vector=c19-1;

    Chauvenete_ddm_mean=mean(Data_input);
    Chauvenete_ddm_std=std(Data_input);

    ddm_clean=Data_input';
end
toc


% xlswrite('excel_file_dWA', height_dWA');
% xlswrite('excel_file_distance', distance_dWA');
% xlswrite('excel_file_ddm', ddm_clean');



% 
% nnn = polyfit(reduced_ZsnrEx,reduced_AmEx_Smth,4)
%  
% fit_is= nnn(1)*(reduced_ZsnrEx).^4+nnn(2)*(reduced_ZsnrEx).^3+nnn(3)*(reduced_ZsnrEx).^2+nnn(4)*(reduced_ZsnrEx) ...
%      +nnn(5)  %%*(reduced_ZsnrEx)+nnn(6);
























