

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%=================== Analytical stuff such as Energy dissipation, Virial, etc. ============================================================%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


PhEx_crop=PhEx(Position_Zc(2,1):Position_Zc(1,1));
DflEx_crop=DflEx(Position_Zc(2,1):Position_Zc(1,1));
ANGLE_CONSERVATIVE_crop=ANGLE_CONSERVATIVE(Position_Zc(2,1):Position_Zc(1,1));
cosine_angle_EX_crop=cosine_angle_EX(Position_Zc(2,1):Position_Zc(1,1));
sine_angle_EX_crop=sine_angle_EX(Position_Zc(2,1):Position_Zc(1,1));


if Simulation==1		%% Not implemented yet
	Edis_T_crop=Edis_T(Position_Zc(2,1):Position_Zc(1,1));
end

%%% rearrange size of shorter vectors
d_min_Incr_d=zeros(length(d_min), 1);   %% the d suffix implies data sorted for analyis with d min sorted 
OMEGA_AM_Dif_d=zeros(length(d_min), 1); 
Fi_Dis_Incr_d=zeros(length(d_min), 1); 
d_min_Incr_d=[d_min_Incr; d_min_Incr(length(d_min)-1)];
OMEGA_AM_Dif_d=[OMEGA_AM_Dif; OMEGA_AM_Dif(length(d_min)-1)];
Fi_Dis_Incr_d=[Fi_Dis_Incr; Fi_Dis_Incr(length(d_min)-1)];


Set_d_min=[d_min, ZsnrEx, AmEx, PhEx, DflEx, ANGLE_CONSERVATIVE, d_min_Incr_d, cosine_angle_EX, sine_angle_EX, OMEGA_AM, OMEGA_AM_Dif_d, Fi_Dis, Fi_Dis_Incr_d, Damping_coefficient];
Set_d_min_crop=[d_min_crop, ZsnrEx_crop, AmEx_crop, PhEx_crop, DflEx_crop, ANGLE_CONSERVATIVE_crop, d_min_Incr_crop, cosine_angle_EX_crop, sine_angle_EX_crop, OMEGA_AM_crop, OMEGA_AM_Dif_crop, Fi_Dis_crop, Fi_Dis_Incr_crop, Damping_coefficient_crop];
			
%%%%%%%%%%%%   Full data sort  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Set_d_min_sorted=sortrows(Set_d_min);
Set_d_min_sorted=flipdim(Set_d_min_sorted,1);   % Sort descending
			
%%%%%%%%%%%%   Croped data sort %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Set_d_min_sorted_Crop=sortrows(Set_d_min_crop);
Set_d_min_sorted_Crop=flipdim(Set_d_min_sorted_Crop,1);   % Sort descending

%%%%%%%%%%%%%%%  Reasign full data for Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d_min_d=Set_d_min_sorted(:,1);	
ZsnrEx_d=Set_d_min_sorted(:,2);		%%% Suffix a stands for analysis: sorted from dmin
AmEx_d=Set_d_min_sorted(:,3);
PhEx_d=Set_d_min_sorted(:,4);
DflEx_d=Set_d_min_sorted(:,5);
ANGLE_CONSERVATIVE_d=Set_d_min_sorted(:,6);
d_min_Incr_d=Set_d_min_sorted(:,7);
cosine_angle_EX_d=Set_d_min_sorted(:,8);
sine_angle_EX_d=Set_d_min_sorted(:,9);
OMEGA_AM_d=Set_d_min_sorted(:,10);
OMEGA_AM_Dif_d=Set_d_min_sorted(:,11);
Fi_Dis_d=Set_d_min_sorted(:,12);
Fi_Dis_Incr_d=Set_d_min_sorted(:,13);
Damping_coefficient_d=Set_d_min_sorted(:,14);

%%%%% Croped  %%%%%%%

d_min_crop_d=Set_d_min_sorted_Crop(:,1);	
ZsnrEx_crop_d=Set_d_min_sorted_Crop(:,2);		%%% Suffix a stands for analysis: sorted from dmin
AmEx_crop_d=Set_d_min_sorted_Crop(:,3);
PhEx_crop_d=Set_d_min_sorted_Crop(:,4);
DflEx_crop_d=Set_d_min_sorted_Crop(:,5);
ANGLE_CONSERVATIVE_crop_d=Set_d_min_sorted_Crop(:,6);
d_min_Incr_crop_d=Set_d_min_sorted_Crop(:,7);
cosine_angle_EX_crop_d=Set_d_min_sorted_Crop(:,8);
sine_angle_EX_crop_d=Set_d_min_sorted_Crop(:,9);
OMEGA_AM_crop_d=Set_d_min_sorted_Crop(:,10);
OMEGA_AM_Dif_crop_d=Set_d_min_sorted_Crop(:,11);
Fi_Dis_crop_d=Set_d_min_sorted_Crop(:,12);
Fi_Dis_Incr_crop_d=Set_d_min_sorted_Crop(:,13);
Damping_coefficient_crop_d=Set_d_min_sorted_Crop(:,14);


	%%%% Sort as function of A ========================================= %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	  %%%% Sort as function of A  ========================================= %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%% Sort as function of A ========================================= %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%%% rearrange size of shorter vectors
d_min_Incr_A=zeros(length(d_min), 1);   %% the d suffix implies data sorted for analyis with d min sorted 
OMEGA_AM_Dif_A=zeros(length(d_min), 1); 
Fi_Dis_Incr_A=zeros(length(d_min), 1); 
d_min_Incr_A=[d_min_Incr; d_min_Incr(length(d_min)-1)];
OMEGA_AM_Dif_A=[OMEGA_AM_Dif; OMEGA_AM_Dif(length(d_min)-1)];
Fi_Dis_Incr_A=[Fi_Dis_Incr; Fi_Dis_Incr(length(d_min)-1)];


%%%%%%

Set_A=[AmEx, ZsnrEx,d_min, PhEx, DflEx, ANGLE_CONSERVATIVE, d_min_Incr_A, cosine_angle_EX, sine_angle_EX, OMEGA_AM, OMEGA_AM_Dif_A, Fi_Dis, Fi_Dis_Incr_A, Damping_coefficient];
Set_A_crop=[AmEx_crop,ZsnrEx_crop, d_min_crop , PhEx_crop, DflEx_crop, ANGLE_CONSERVATIVE_crop, d_min_Incr_crop, cosine_angle_EX_crop, sine_angle_EX_crop, OMEGA_AM_crop, OMEGA_AM_Dif_crop, Fi_Dis_crop, Fi_Dis_Incr_crop, Damping_coefficient_crop];
			
%%%%%%%%%%%%   Full data sort  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Set_A_sorted=sortrows(Set_A);
Set_A_sorted=flipdim(Set_A_sorted,1);   % Sort descending
			
%%%%%%%%%%%%   Croped data sort %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Set_A_sorted_Crop=sortrows(Set_A_crop);
Set_A_sorted_Crop=flipdim(Set_A_sorted_Crop,1);   % Sort descending

%%%%%%%%%%%%%%%  Reasign full data for Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AmEx_A=Set_A_sorted(:,1);	
ZsnrEx_A=Set_A_sorted(:,2);		%%% Suffix a stands for analysis: sorted from A
d_min_A=Set_A_sorted(:,3);
PhEx_A=Set_A_sorted(:,4);
DflEx_A=Set_A_sorted(:,5);
ANGLE_CONSERVATIVE_A=Set_A_sorted(:,6);
d_min_Incr_A=Set_A_sorted(:,7);
cosine_angle_EX_A=Set_A_sorted(:,8);
sine_angle_EX_A=Set_A_sorted(:,9);
OMEGA_AM_A=Set_A_sorted(:,10);
OMEGA_AM_Dif_A=Set_A_sorted(:,11);
Fi_Dis_A=Set_A_sorted(:,12);
Fi_Dis_Incr_A=Set_A_sorted(:,13);
Damping_coefficient_A=Set_A_sorted(:,14);

%%%%% Croped  %%%%%%%

AmEx_crop_A=Set_A_sorted_Crop(:,1);	
ZsnrEx_crop_A=Set_A_sorted_Crop(:,2);		%%% Suffix a stands for analysis: sorted from A
d_min_crop_A=Set_A_sorted_Crop(:,3);
PhEx_crop_A=Set_A_sorted_Crop(:,4);
DflEx_crop_A=Set_A_sorted_Crop(:,5);
ANGLE_CONSERVATIVE_crop_A=Set_A_sorted_Crop(:,6);
d_min_Incr_crop_A=Set_A_sorted_Crop(:,7);
cosine_angle_EX_crop_A=Set_A_sorted_Crop(:,8);
sine_angle_EX_crop_A=Set_A_sorted_Crop(:,9);
OMEGA_AM_crop_A=Set_A_sorted_Crop(:,10);
OMEGA_AM_Dif_crop_A=Set_A_sorted_Crop(:,11);
Fi_Dis_crop_A=Set_A_sorted_Crop(:,12);
Fi_Dis_Incr_crop_A=Set_A_sorted_Crop(:,13);
Damping_coefficient_crop_A=Set_A_sorted_Crop(:,14);



			