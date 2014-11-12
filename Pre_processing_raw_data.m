%% ======================COLLECT curves =================================================================================== %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
	Collect_curves_file;   % collects curves Experimental and simulation 

	%%%%%%%%%========= END COLLECT CURVES   ====================================================================================%%%%%%%%%%%%%%%%%%%%%%%%%


	%%==============Sorting zc cantilever separation from large to small not changing order of elements, i.e. raw data===========%%%%%%%%%%%%%%%%%%%
  
    Sorting_zc_cantilever_separation;    % this will force first element zc in vector to be zc>>A0 and last Zc<<A0 no change in data
	
	%%------------------ END CANTILEVER SEPARATION ORDER, i.e. raw data ===========================================================%%%%%%%%%%%%%%%%%%%
    
	%%%%%  Amax and deflection in nm ==========================================================================================%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
	% Experimental only since simulation is in nm
	AmEx=AmEx*AmpInvOLS*1e-9;
	DflEx=DflEx*AmpInvOLS*1e-9/1.09;
	
	%%% IMPORTANT, from here on Amplitude and deflection are in nm ====================================================%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	%%%%% ================ Ploting dmin =======================================================================%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if plot_raw_curves==1
		figure (count_figures+86) 
		plot(ZsnrEx, ZsnrEx-AmEx, '.k', 'Markersize',3, 'displayname','D_min');    %% This is dmin
		title(' Dmin versus zc: raw','fontsize',12)
		xlabel('Zc','fontsize',14) 
		ylabel('dmin','fontsize',14)
	end
	%%%%% ================ END Ploting dmin RAW =======================================================================%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Preparing DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
	%%%% remove first and last part of vectors ============================================================================================ %%%%%%%%%%%%%%%
	
    Set_remove=[ZsnrEx, AmEx, PhEx, DflEx];
    Set_remove=Set_remove(remove_start:end-remove_end,:);   %% Because begining and end might induce errors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ZsnrEx=Set_remove(:,1);
    AmEx=Set_remove(:,2);
    PhEx=Set_remove(:,3);
    DflEx=Set_remove(:,4);
  
    
    
	
	%%%%%% Allow to crop curve %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	if crop_initial_data==1
		
        figure (90)
        plot(ZsnrEx,AmEx, '.k', 'Markersize',M_size, 'displayname','Amplitud A'); 
        hold on
        [x1,y1]=ginput(2);

        Pck_ini=[x1 , y1];
        Pck_ini=sort(Pck_ini,1); 


        for ii=1:1:2

            Zc_values_ini(ii,1)=Pck_ini(ii,1); 
            Minimum_pick_ini(ii,1)=Pck_ini(ii,2);
            plot(Zc_values_ini(ii,1),  Minimum_pick_ini(ii,1),'Vk', 'Markersize',5)    
            Minimum_error_dumb_ini=MultiplierError*min(abs(ZsnrEx-Zc_values_ini(ii,1)));
            Minimum_error_ini(ii,1)=Minimum_error_dumb_ini;

            Position_Zc_dumb_ini=find(abs(ZsnrEx-Zc_values_ini(ii,1))<=abs(Minimum_error_dumb_ini),1, 'last');
            Position_Zc_ini(ii,1)=Position_Zc_dumb_ini; % this is the position of the vector element
            Zc_found_ini(ii,1)=ZsnrEx(Position_Zc_dumb_ini);
            plot(Zc_found_ini(ii,1),  AmEx(Position_Zc_dumb_ini),'Vk', 'Markersize',5)
        end  

        if Simulation==0

            Set_to_crop=[ZsnrEx, AmEx, PhEx, DflEx];
            Set_to_crop=Set_to_crop(Position_Zc_ini(2,1):Position_Zc_ini(1,1),:);

            ZsnrEx=Set_to_crop(:,1);
            AmEx=Set_to_crop(:,2);
            PhEx=Set_to_crop(:,3);
            DflEx=Set_to_crop(:,4);

        else
            Set_to_crop=[ZsnrEx, AmEx, PhEx, DflEx, Edis_T];
            Set_to_crop=Set_to_crop(Position_Zc_ini(1,1):Position_Zc_ini(2,1),:);

            ZsnrEx=Set_to_crop(:,1);
            AmEx=Set_to_crop(:,2);
            PhEx=Set_to_crop(:,3);
            DflEx=Set_to_crop(:,4);
            Edis_T=Set_to_crop(:,5);
        end
	
	end	

	
	%%%% End remove first and last part of vectors ============================================================================= %%%%%%%%%%%%%%%
	

	
	%%%% ============ Sorting zc values (OPTIONAL) from larger to smaller ==========================================================%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
    if Sorting_zc_larger_to_smaller==1    %% This activated might induce errors
        Sorting_zc_larger_to_smaller_file;
    end
    
    if Removing_repeat_zc==1  %  might induce errors due to lack of data points, in particular if they are not sorted from large to small zc
		Removing_repeat_zc_file;
            
    end  %% end of removing same Zc-------------------------------------------------------------------------------------------------------%%
        
	%%%%%%%%%% end sorting zc and/or removing zc ================================================================================== %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    