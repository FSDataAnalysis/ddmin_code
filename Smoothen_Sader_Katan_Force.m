%%%%%%%%%%% REMOVING ZEROES; INFS AND COMPLEX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%% Removing NANS %%%%%%%%%%%%%%%%%%%%%%%%%%

Position_FTS_K_S_NAN=find(isnan(Fts_cons_crop_raw(:)));  %%% removing NANS

Number_Nans_K_S=length(Position_FTS_K_S_NAN);			 




Fts_cons_crop_Complex=Fts_cons_crop_raw;
d_min_crop_Complex=d_min_crop;

Fts_cons_crop_Complex(Position_FTS_K_S_NAN)=[];
d_min_crop_Complex(Position_FTS_K_S_NAN)=[];


%%%%%%%%%%%% Removing INFS %%%%%%%%%%%%%%%%%%%%%%%%%%
														 
														
Position_FTS_K_S_inf=find(real(Fts_cons_crop_Complex)==Inf);  % %%% removing infs

Number_infs_K_S=length(Position_FTS_K_S_inf);			 

Fts_cons_crop_Complex(Position_FTS_K_S_inf)=[];
d_min_crop_Complex(Position_FTS_K_S_inf)=[];

Position_FTS_K_S_minus_inf=find(real(Fts_cons_crop_Complex)==-Inf);  % %%% removing infs

Number_minus_infs_K_S=length(Position_FTS_K_S_minus_inf);			 


Fts_cons_crop_Complex(Position_FTS_K_S_minus_inf)=[];
d_min_crop_Complex(Position_FTS_K_S_minus_inf)=[];

%%%%%%%%%%%% Removing zeros %%%%%%%%%%%%%%%%%%%%%%%%%%

Position_FTS_K_S_zeros=find(real(Fts_cons_crop_Complex)==0);   % %%% removing zeros

Number_zeros_K_S=length(Position_FTS_K_S_zeros);			 


Fts_cons_crop_Complex(Position_FTS_K_S_zeros)=[];
d_min_crop_Complex(Position_FTS_K_S_zeros)=[];

%%%%%%%%%%%  Removing complex %%%%%%%%%%%%%%%%%%%%%%%%%
														
														
Position_FTS_Complex=find(imag(Fts_cons_crop_Complex)~=0);

Number_Complex_K_S=length(Position_FTS_Complex);

Fts_cons_K_T_Real=Fts_cons_crop_Complex;
d_min_K_T_Real=d_min_crop_Complex;
Fts_cons_K_T_Real(Position_FTS_Complex)=[];
d_min_K_T_Real(Position_FTS_Complex)=[];

%%%%%%%%%%% FINAL Fts with real part of complex numbers only  %%%%%%%%%%%%%%%%%%%%%%%%%

Fts_cons_K_T_Final=real(Fts_cons_crop_Complex);
d_min_K_T_Final=real(d_min_crop_Complex);

%%%%%%%%%%% Plotting figure for conservative and effective coefficients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure (count_figures+11)
hold on
%plot(d_min_K_T_Final,Fts_cons_K_T_Final,'.r','Markersize',M_size,'displayname','Conservative force: Sader-Katan:raw'); NOT ZEROED



%%%%%%%%%%% PLOTING SMOOTH DATA ETC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		S_K_set=[d_min_K_T_Final, Fts_cons_K_T_Final];
		S_K_set_sorted=sortrows(S_K_set);
 		S_K_set_sorted=flipdim(S_K_set_sorted,1);   % Sort descending
		d_min_crop_sorted=S_K_set_sorted(:,1);
		Fts_cons_crop_sorted=S_K_set_sorted(:,2);
		
		Fts_cons_crop_sorted_ss=smooth(Fts_cons_crop_sorted,ss_Fts_cons_S_K,'rloess');  %ss_Fts_cons_S_K
		%%%%%%%%%---------------------------- get dmin to zero at Fad -------------------------------------------------------------------%%%%%%%%%%%%%%

		F_adhesion=min(Fts_cons_crop_sorted_ss);
		F_maximum=max(Fts_cons_crop_sorted_ss);    
	
        Minimum_error_Fa=MultiplierError*min(abs(Fts_cons_crop_sorted_ss-F_adhesion));        
        Position_Fa=find(abs(Fts_cons_crop_sorted_ss-F_adhesion)<=abs(Minimum_error_Fa),1, 'last');
        D_min_zero_found=d_min_crop_sorted(Position_Fa);
		
		d_min_crop_sorted_zeroed=d_min_crop_sorted-d_min_crop_sorted(Position_Fa);
		d_min_crop_zeroed=d_min_crop-d_min_crop_sorted(Position_Fa);
		
		
		%%%%% Plot zeroed DMIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		figure (count_figures+11)    
		
		hold on
		plot(d_min_crop_zeroed,Fts_cons_crop_raw,'.r','Markersize',M_size,'displayname','Conservative force: Sader-Katan:raw'); 
		hold on
		plot(d_min_crop_sorted_zeroed,Fts_cons_crop_sorted_ss,'.b','Markersize',M_size, 'displayname','Conservative force: Sader-Katan: Smooth');
		box on
		title('Force Sader and damping coefficient versus zc','fontsize',12)
		xlabel('dmin','fontsize',14) 
		ylabel('Conservative force and damping','fontsize',14); 
	
				
		axis([(d_min_crop_sorted_zeroed(end)-1e-9) d_min_crop_sorted_zeroed(1) -abs(2*F_adhesion) 4*(abs(F_adhesion))])
		
		
		%%%% Smoothen Damping coefficient %%%%%%%%%%%%%%%%%%%%%
		
		hold on
		plot(d_min_crop_zeroed,Damping_coefficient_crop,'.k','Markersize',M_size, 'displayname','Damping coefficient raw');
		
		Set_Damp_Sorted=[d_min_crop_zeroed, Damping_coefficient_crop];
		Set_Damp_Sorted=sortrows(Set_Damp_Sorted);
		Set_Damp_Sorted=flipdim(Set_Damp_Sorted,1); 
		d_min_Damp=Set_Damp_Sorted(:,1);
		Damping_coefficient_Sorted=Set_Damp_Sorted(:,2);
		Damping_coefficient_crop_ss=smooth(Damping_coefficient_Sorted,ss_Damping_Coeff,'rloess');						% What ends us with S is smooth from final curve not smoothed data 
		plot(d_min_Damp,Damping_coefficient_crop_ss,'.c','Markersize',M_size, 'displayname','Damping coefficient:Smooth');	
	
		Damping_coefficient_crop_ss_max=max(Damping_coefficient_crop_ss);
		
		text(0,1*abs(F_adhesion), ['Max damping :' num2str(Damping_coefficient_crop_ss_max), ' Ns/m'],'fontsize',12)
		text(0,0.6*abs(F_adhesion), ['F adhesion:' num2str(F_adhesion), ' N'],'fontsize',12)

		
%%%%% End fig 11 %%%%%%%%%%
%%%%%%%%%%%%%%%%%% INTEGRAL PART %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


d_min_crop_dis=d_min_crop(1:end-1);
Fts_dis_Complex=Fts_dis_crop_raw;

%%%%%%%NANS 
Fts_dis_Complex_NAN=find(isnan(Fts_dis_Complex(:)));  %%% removing NANS
Fts_dis_Complex(Fts_dis_Complex_NAN)=[];
d_min_crop_dis(Fts_dis_Complex_NAN)=[];

%%%%INFS%%%%%
Position_dis_inf=find(real(Fts_dis_Complex)==Inf);  % %%% removing infs
Fts_dis_Complex(Position_dis_inf)=[];
d_min_crop_dis(Position_dis_inf)=[];

Position_dis_inf_minus=find(real(Fts_dis_Complex)==-Inf);  % %%% removing infs
Fts_dis_Complex(Position_dis_inf_minus)=[];
d_min_crop_dis(Position_dis_inf_minus)=[];

%%% ZEROS %%%

Position_dis_zeros=find(real(Fts_dis_Complex)==0);   % %%% removing zeros
Fts_dis_Complex(Position_dis_zeros)=[];
d_min_crop_dis(Position_dis_zeros)=[];

%%%%%%  FINAL REAL NUMBERS %%%%
Fts_K_S_dis=real(Fts_dis_Complex);
d_min_K_S_dis=real(d_min_crop_dis);


%%%%%%  SORTING FROM SMALLER TO LARGER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


S_K_set_dis=[d_min_K_S_dis, Fts_K_S_dis];
S_K_set_dis_sorted=sortrows(S_K_set_dis);
S_K_set_dis_sorted=flipdim(S_K_set_dis_sorted,1);   % Sort descending
d_min_sorted_dis=S_K_set_dis_sorted(:,1);
Fts_sorted_dis=S_K_set_dis_sorted(:,2);

%%%%%%%%%%%%%%%%%%%%%%%%% PLOTING CONSERVATIVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure (count_figures+12)
hold on

plot(d_min_crop_zeroed,Fts_cons_crop_raw,'.r','Markersize',M_size,'displayname','Conservative force: Sader-Katan:raw'); 
hold on
%Fts_dis_crop_raw_s=smooth(Fts_dis_crop_raw,s_Fts_dis_s,'rloess'); %	
%plot(d_min_crop(1:end-1)-D_min_zero_found,Fts_dis_crop_raw_s,'.m','Markersize',M_size,'displayname','Dissipative coefficient: Integral_smooth');
plot(d_min_crop_sorted_zeroed,Fts_cons_crop_sorted_ss,'.b','Markersize',M_size, 'displayname','Conservative force: Sader-Katan: Smooth');
box on
title('Force Sader and Integral coefficient versus zc','fontsize',12)
xlabel('dmin','fontsize',14) 
ylabel('Conservative force and Integral damping','fontsize',14); 
				
axis([(d_min_crop_sorted_zeroed(end)-1e-9) d_min_crop_sorted_zeroed(1) -abs(2*F_adhesion) 1.2*F_maximum])


%%%%%%   INTEGRAL DISSIPATIVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(d_min_crop(1:end-1)-D_min_zero_found,Fts_dis_crop_raw,'.k','Markersize',M_size, 'displayname','Dissipative coefficient: Sader Integral:raw');
Fts_sorted_dis_ss=smooth(Fts_sorted_dis, s_F_dis,'rlowess');  % ss_Fts_cons_S_K
plot(d_min_sorted_dis-D_min_zero_found,Fts_sorted_dis_ss,'.c','Markersize',M_size, 'displayname','Dissipative coefficient: Sader Integral:Smooth');





