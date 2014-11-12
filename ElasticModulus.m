% 
% Delta_ind=4.4e-9;
% F_rep=7e-9;
% F_a=1.5e-9;
% NetForce=F_rep+F_a;
% 
% E_effective=3*NetForce/(4*(R^(0.5))*(Delta_ind^(1.5)))
% clear all;
% close all;
% load exampleElastic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARAMETERS %%%%%

Select_adhesion_point_manually_ElasticM=1;
NonLinear_ElasticM=1;  %% Calculates Elastic Modulus as a fit
CONFIDENCE_INTERVAL_ElasticM=0.05;
M_size_E=7;
%%%%%%%%%%%%%%% Clean figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



figure (count_figures+101)
hold on
box on


if Sader_Katan==1
		
		hold on
		plot(D_M_ZEROED_CROP,Fts_S_K_cons_norm,'.r','Markersize',M_size,'displayname','Conservative force: Sader-Katan:raw'); 
		plot(d_min_crop_sorted_zeroed,Fts_S_K_cons_norm_ss,'.b','Markersize',M_size, 'displayname','Conservative force: Sader-Katan: Smooth');
		axis([ min_d_plot  max_d_plot -1.2 2])	
		text(1e-9,0.8, ['A0 :' num2str(A0*1e9) '  nm'],'fontsize',12);
		text(1e-9,0.4, ['F adhesion:' num2str(F_adhesion), ' N'],'fontsize',12)
else

		text(1e-9,0.6, ['A0 :' num2str(A0*1e9) '  nm'],'fontsize',12);
		axis([ d_min(1) d_min(end)  0 1.2])
end



%% pick points to crop %%%%%%

%%%% e stands for elastic modulus 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Select_adhesion_point_manually_ElasticM==1   % This allows you pick three points, i.e. adhesion and two for range
    title('Elastic Modulus: Pick Range Indentation and zero load','fontsize',12)
    %This section of the code has been modified to pick the points
    %automatically. Intersectline1 selects is the leftmost part of the range at
    %a specific y value
%     intersectline1=[-3e-009 0.5;1e-009 0.5];   %intersectline1=[-5e-009 1.5;5e-009 1.5];
%     intercepts1=polyxpoly(intersectline1(:,1),intersectline1(:,2),d_min_crop_sorted_zeroed,Fts_S_K_cons_norm_ss);
%     f1=length(intercepts1);
%     %we choose the rightmost intercept
%     r1=intercepts1(f1,1);
    %The rightmost part of the range coincides with zero force.
%     intersectline2=[-1e-009 0;1e-009 0];
%     intercepts2=polyxpoly(intersectline2(:,1),intersectline2(:,2),d_min_crop_sorted_zeroed,Fts_S_K_cons_norm_ss);
%     f2=length(intercepts2);
% %     r2=intercepts2(f2,1);
%     r2=intercepts2(f2,1);

    dumb_pos=Position_Fa;
    Pck_e=[d_min_crop_sorted_zeroed(end-subtract_EModulus) Fts_S_K_cons_norm_ss(end-subtract_EModulus); d_min_crop_sorted_zeroed(dumb_pos+add_EModulus) Fts_S_K_cons_norm_ss(dumb_pos+add_EModulus); d_min_crop_sorted_zeroed(dumb_pos) Fts_S_K_cons_norm_ss(dumb_pos)];

    %To select points manually, use the following code instead
    %[q1,q2, q3]=ginput(3);

    %Pck_e=[q1 , q2, q3];
    %Pck_e=sort(Pck_e,1); 

    this=0;

    %%%%%%%%%%%%%%%%%%%%%    Find points %%%%%%%%%%%%%%
    for ii=1:1:3

            d_values_e(ii,1)=Pck_e(ii,1); 
            Minimum_pick_e(ii,1)=Pck_e(ii,2);
            plot(d_values_e(ii,1),  Minimum_pick_e(ii,1),'Vk', 'Linewidth', 2, 'Markersize',6);

            Minimum_error_dumb_e=MultiplierError*min(abs(d_min_crop_sorted_zeroed-d_values_e(ii,1)));
            Minimum_error_e(ii,1)=Minimum_error_dumb_e;
            Position_d_dumb_e=find(abs(d_min_crop_sorted_zeroed-d_values_e(ii,1))<=abs(Minimum_error_dumb_e),1, 'first');
            Position_d_dumb_ElasticM(ii,1)=Position_d_dumb_e;

            Energy_found_e(ii,1)=Fts_S_K_cons_norm_ss(Position_d_dumb_e);
            d_found_e(ii,1)=d_min_crop_sorted_zeroed(Position_d_dumb_e);
            plot(d_found_e(ii,1),  Fts_S_K_cons_norm_ss(Position_d_dumb_e),'Vk','Linewidth', 4, 'Markersize',6)

    end  


    % %%%%%% Prepare vectors to calculate K %%%%%%%%%%%%%%%%%%%%
    EE1=Position_d_dumb_ElasticM(2,1);
    EE2=Position_d_dumb_ElasticM(1,1); 

    Position_Adhesion=Position_d_dumb_ElasticM(3,1);
    Picked_Adhesion=Fts_S_K_cons_norm_ss(Position_Adhesion);
    Picked_Adhesion_abs=abs(Picked_Adhesion*F_adhesion);

    e_Force=Fts_S_K_cons_norm_ss;
    e_Indentation=d_min_crop_sorted_zeroed;

    % %%% Standard  International Units %%%%%%

    e_Force_IU=Fts_S_K_cons_norm_ss*(abs(F_adhesion));

    e_Indentation_IU=abs(d_min_crop_sorted_zeroed(EE1:EE2)-d_min_crop_sorted_zeroed(Position_Adhesion));

    e_load_IU=e_Force_IU(EE1:EE2)+abs(F_adhesion);

    Constant_e=3/4*R^(-0.5);


    Elastic_Modulus_Analytic=Constant_e*(e_load_IU.*(e_Indentation_IU.^(-1.5)));
    Elastic_Modulus_Analytic_dumb=Elastic_Modulus_Analytic;

    % %%%% REMOVE NANS ETC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%  curve raw %%%%
        DUMMY_YY=zeros(length(Elastic_Modulus_Analytic));
        DUMMY_XX=zeros(length(e_Indentation_IU));
        DUMMY_YY=Elastic_Modulus_Analytic;
        DUMMY_XX=e_Indentation_IU;

        GetRidOfNANSINFS; %%% This gets rid of infs, NANS,complex etc, and returns real and numbers only as DUMMY_YY_real and DUMMY_XX_real
        E_real=DUMMY_YY_real;
        Indentation_real=DUMMY_XX_real;



    Elastic_Modulus_Analytic_mean=mean(E_real);
    Elastic_Modulus_Analytic_std=std(E_real);

    %%%% remove outliers %%%%%%%

        Data_input=E_real;

        Chauvenete;
        Outliers_E_vector=c19;

    Chauvenete_ElasticM_mean=mean(Data_input);
    Chauvenete_ElasticM_std=std(Data_input);



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if NonLinear_ElasticM==1 %%% Non linear fit for elastic modulus


    y=e_load_IU;
    dy = y*CONFIDENCE_INTERVAL_ElasticM;

    xx_e=4/3*(R^(0.5))*(e_Indentation_IU.^(1.5));
    x=[xx_e];    % 

    %formula converted to
    %The inline version
    func = inline('p(1)*x(:,1)','p','x');   % %%% Here p is the elastic modulus

    
    
    p0 = Chauvenete_ElasticM_mean;

    %To detect the sensitivity of the fit to starting parameter guess,
    %the fit is run a number of times.
    %each fit is plotted and each parameter plotted as a histogram
    Nrepeat=20;
    %each parameter is varied by a normal distribution with
    %mean equal to the starting guess and std.dev. equal to
    %sd*mean
    sd = 0.3;
    %histogram zoom factor (how many std dev to show)
    zfactor = 2;
    %parameter outlier cuttoff: lowest and highest N estimates are removed
    outcut=10;
    %========================================================
    %END settable inputs
    %========================================================
    terms=1;
    %list of all parameter outputs to use in histogram
    pList=zeros(Nrepeat,terms);

    for rep =1:Nrepeat
        %     rep

        %form the new randomized start vector
        p = [p0(1)*(1+sd*randn)];
        %do the fit
        [p,r,j] = nlinfit(x,y,func,p);
        %copy fit to list
        pList(rep,:) = p';

        %get parameter errors
        c95 = nlparci(p,r,j);
        [yp, ci] = nlpredci(func,x,p,r,j);

        %plot the fit 1 parameter only
        %if Plot_process_quantification==1 
            figure (count_figures+102)
            hold on 
            box on
            title(' Elastic Modulus: NONLINEAR FIT','fontsize',12)
            xlabel('Indentation','fontsize',14) 
            ylabel('Force: NonLinear FIT Netwons','fontsize',14)

            errorbar(e_Indentation_IU,func(p,x),ci,ci,'b-');			% fit

            hold on

            errorbar(e_Indentation_IU,y,dy,dy,'ro')		%% real 


    end

    figure (count_figures+103)
    hold on 
    box on

    output_ElasticM_mean=mean(pList);
    output_ElasticM_std=std(pList);
    format shortEng

    figure (count_figures+103)
    hold on 
    box on

    Force_nl_reconstructed=output_ElasticM_mean*xx_e;

    Force_nl_reconstructed_net=Force_nl_reconstructed+F_adhesion;
    Force_nl_reconstructed_net_norm=Force_nl_reconstructed_net/Picked_Adhesion_abs;

    plot(-e_Indentation_IU,Force_nl_reconstructed_net_norm,'Vk','Markersize',M_size_E, 'displayname',' NONLINEAR FIT From ElasticM (DMT): Force');

    text(1e-9,-1.4, ['ElasticM (DMT) NN:' num2str(output_ElasticM_mean*1e-9) '  GPa'],'fontsize',12);
    text(1e-9,-1.8, ['ElasticM (DMT) NN (std) :' num2str(output_ElasticM_std*1e-9) '  GPa'],'fontsize',12);
    axis([ min_d_plot  max_d_plot -2 2])

    end   % non linear fit

    figure (count_figures+103)
    hold on 
    box on
    DD_subtract=d_min_crop_sorted_zeroed(Position_Adhesion);
    plot(D_M_ZEROED_CROP-DD_subtract,Fts_S_K_cons_norm*abs(F_adhesion)/Picked_Adhesion_abs,'.r','Markersize',M_size,'displayname','Conservative force: Sader-Katan:raw'); 
    plot(d_min_crop_sorted_zeroed-DD_subtract,Fts_S_K_cons_norm_ss*abs(F_adhesion)/Picked_Adhesion_abs,'.b','Markersize',M_size, 'displayname','Conservative force: Sader-Katan: Smooth');

    title(' Elastic Modulus: FIT','fontsize',12)
    axis([ min_d_plot  max_d_plot -2 2.2])	

    text(1e-9,1.2, ['Range Indentation :' num2str(e_Indentation_IU(1)*1e9) ' to ' num2str(e_Indentation_IU(end)*1e9) '  nm'],'fontsize',12);
    text(1e-9,0.8, ['A0 :' num2str(A0*1e9) '  nm'],'fontsize',12);
    text(1e-9,0.4, ['F adhesion:' num2str(F_adhesion), ' N'],'fontsize',12)

    %%%%%% Plot Fit from Analytic 

    Force_Analytic_FIT=4/3*Chauvenete_ElasticM_mean*(R^(0.5))*(e_Indentation_IU.^(1.5));
%     plot(-e_Indentation_IU,Force_Analytic_FIT/Picked_Adhesion_abs +F_adhesion/Picked_Adhesion_abs,'om','Markersize',M_size_E, 'displayname',' Analytic FIT From ElasticM (DMT): FORCE');
% 
%     text(1e-9,-0.4, ['ElasticM (DMT) AL:' num2str(Chauvenete_ElasticM_mean*1e-9) '  GPa'],'fontsize',12);
%     text(1e-9,-0.8, ['ElasticM (DMT) AL (std) :' num2str(Chauvenete_ElasticM_std*1e-9) '  GPa'],'fontsize',12);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%  Calculating error Elastic modulus anlytic  fit and non linear  %%%%%%%%%


    figure (count_figures+103)
    hold on 
    box on
    Force_Analytic_FIT_total_norm=Force_Analytic_FIT/Picked_Adhesion_abs +F_adhesion/Picked_Adhesion_abs;
    Force_total_norm =e_load_IU/Picked_Adhesion_abs+F_adhesion/Picked_Adhesion_abs;
    % Difference_squared_force_fit=(((Force_total_norm-Force_Analytic_FIT_total_norm)./Force_total_norm).^2);
    Difference_squared_force_fit=(((Force_total_norm-Force_Analytic_FIT_total_norm)/max(abs(Force_total_norm))).^2);
% 
    Data_input=Difference_squared_force_fit;

    Chauvenete;
    Sum_squared_force_fit=sum(Data_input);
    force_fit_points=length(Data_input);

    Standard_deviation_error_force_fit=floor(sqrt(1/force_fit_points*(Sum_squared_force_fit))*100);
% 
%     text(0,1.6, ['ERROR ANALTIC FIT (STD): ' num2str(Standard_deviation_error_force_fit), ' %'],'fontsize',15);


    if NonLinear_ElasticM==1

    % 	Difference_squared_force_nl_fit=(((Force_total_norm-Force_nl_reconstructed_net_norm)./Force_total_norm).^2);
        Difference_squared_force_nl_fit=(((Force_total_norm-Force_nl_reconstructed_net_norm)/max(abs(Force_total_norm))).^2);


        Data_input=Difference_squared_force_nl_fit;

        Chauvenete;
        Sum_squared_force_nl_fit=sum(Data_input);
        force_nl_fit_points=length(Data_input);

        Standard_deviation_error_force_nl_fit=floor(sqrt(1/force_nl_fit_points*(Sum_squared_force_nl_fit))*100);

        text(0,2, ['ERROR NL FIT (STD): ' num2str(Standard_deviation_error_force_nl_fit), ' %'],'fontsize',15);


    end




    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    else  %%% Takes automatic adhesion point (pick range only)
    title('Elastic Modulus: Pick Range Indentation','fontsize',12)	

    [q1,q2]=ginput(2);

    Pck_e=[q1 , q2];
    Pck_e=sort(Pck_e,1); 



    %%%%%%%%%%%%%%%%%%%%%    Find points %%%%%%%%%%%%%%
    for ii=1:1:2

            d_values_e(ii,1)=Pck_e(ii,1); 
            Minimum_pick_e(ii,1)=Pck_e(ii,2);
            plot(d_values_e(ii,1),  Minimum_pick_e(ii,1),'Vk', 'Linewidth', 2, 'Markersize',6);

            Minimum_error_dumb_e=MultiplierError*min(abs(d_min_crop_sorted_zeroed-d_values_e(ii,1)));
            Minimum_error_e(ii,1)=Minimum_error_dumb_e;
            Position_d_dumb_e=find(abs(d_min_crop_sorted_zeroed-d_values_e(ii,1))<=abs(Minimum_error_dumb_e),1, 'first');
            Position_d_dumb_ElasticM(ii,1)=Position_d_dumb_e;

            force_found_e(ii,1)=Fts_S_K_cons_norm_ss(Position_d_dumb_e);
            d_found_e(ii,1)=d_min_crop_sorted_zeroed(Position_d_dumb_e);
            plot(d_found_e(ii,1),  Fts_S_K_cons_norm_ss(Position_d_dumb_e),'Vk','Linewidth', 4, 'Markersize',6)

    end  

    % %%%%%% Prepare vectors to calculate K %%%%%%%%%%%%%%%%%%%%
    EE1=Position_d_dumb_ElasticM(2,1);
    EE2=Position_d_dumb_ElasticM(1,1); 

    e_Force=Fts_S_K_cons_norm_ss;
    e_Indentation=d_min_crop_sorted_zeroed;

    % %%% Standard  International Units %%%%%%

    e_Force_IU=Fts_S_K_cons_norm_ss*(abs(F_adhesion));

    e_Indentation_IU=abs(d_min_crop_sorted_zeroed(EE1:EE2));

    e_load_IU=e_Force_IU(EE1:EE2)+abs(F_adhesion);

    Constant_e=3/4*R^(-0.5);


    Elastic_Modulus_Analytic=Constant_e*(e_load_IU.*(e_Indentation_IU.^(-1.5)));
    Elastic_Modulus_Analytic_dumb=Elastic_Modulus_Analytic;

    % %%%% REMOVE NANS ETC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%  curve raw %%%%
        DUMMY_YY=zeros(length(Elastic_Modulus_Analytic));
        DUMMY_XX=zeros(length(e_Indentation_IU));
        DUMMY_YY=Elastic_Modulus_Analytic;
        DUMMY_XX=e_Indentation_IU;

        GetRidOfNANSINFS; %%% This gets rid of infs, NANS,complex etc, and returns real and numbers only as DUMMY_YY_real and DUMMY_XX_real
        E_real=DUMMY_YY_real;
        Indentation_real=DUMMY_XX_real;



    Elastic_Modulus_Analytic_mean=mean(E_real);
    Elastic_Modulus_Analytic_std=std(E_real);

    %%%%% remove outliers %%%%%%%

        Data_input=E_real;

        Chauvenete;
        Outliers_E_vector=c19;

    Chauvenete_ElasticM_mean=mean(Data_input);
    Chauvenete_ElasticM_std=std(Data_input);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if NonLinear_ElasticM==1 %%% Non linear fit for elastic modulus


    y=e_load_IU;
    dy = y*CONFIDENCE_INTERVAL_ElasticM;

    xx_e=4/3*(R^(0.5))*(e_Indentation_IU.^(1.5));
    x=[xx_e];    % 

    %formula converted to
    %The inline version
    func = inline('p(1)*x(:,1)','p','x');   % %%% Here p is the elastic modulus

    p0 = [Chauvenete_ElasticM_mean];

    %To detect the sensitivity of the fit to starting parameter guess,
    %the fit is run a number of times.
    %each fit is plotted and each parameter plotted as a histogram
    Nrepeat=20;
    %each parameter is varied by a normal distribution with
    %mean equal to the starting guess and std.dev. equal to
    %sd*mean
    sd = 0.3;
    %histogram zoom factor (how many std dev to show)
    zfactor = 2;
    %parameter outlier cuttoff: lowest and highest N estimates are removed
    outcut=10;
    %========================================================
    %END settable inputs
    %========================================================
    terms=1;
    %list of all parameter outputs to use in histogram
    pList=zeros(Nrepeat,terms);

    for rep =1:Nrepeat
        %     rep

        %form the new randomized start vector
        p = [p0(1)*(1+sd*randn)];
        %do the fit
        [p,r,j] = nlinfit(x,y,func,p);
        %copy fit to list
        pList(rep,:) = p';

        %get parameter errors
        c95 = nlparci(p,r,j);
        [yp, ci] = nlpredci(func,x,p,r,j);

        %plot the fit 1 parameter only
        %if Plot_process_quantification==1 
            figure (count_figures+102)
            hold on 
            box on
            title(' Elastic Modulus: NONLINEAR FIT','fontsize',12)
            xlabel('Indentation','fontsize',14) 
            ylabel('Force: NonLinear FIT Netwons','fontsize',14)

            errorbar(e_Indentation_IU,func(p,x),ci,ci,'b-');

            hold on

            errorbar(e_Indentation_IU,y,dy,dy,'ro')


    end

    figure (count_figures+103)
    hold on 
    box on

    output_ElasticM_mean=mean(pList);
    output_ElasticM_std=std(pList);
    format shortEng

    figure (count_figures+103)
    hold on 
    box on

    Force_nl_reconstructed=output_ElasticM_mean*xx_e;


    plot(-e_Indentation_IU,Force_nl_reconstructed/abs(F_adhesion)-1,'Vk','Markersize',M_size_E, 'displayname',' NONLINEAR FIT From ElasticM (DMT): Force');

    text(1e-9,-1, ['ElasticM (DMT) NN:' num2str(output_ElasticM_mean*1e-9) '  GPa'],'fontsize',12);
    text(1e-9,-0.8, ['ElasticM (DMT) NN (std) :' num2str(output_ElasticM_std*1e-9) '  GPa'],'fontsize',12);

    end

    figure (count_figures+103)
    hold on 
    box on
    plot(D_M_ZEROED_CROP,Fts_S_K_cons_norm,'.r','Markersize',M_size,'displayname','Conservative force: Sader-Katan:raw'); 
    plot(d_min_crop_sorted_zeroed,Fts_S_K_cons_norm_ss,'.b','Markersize',M_size, 'displayname','Conservative force: Sader-Katan: Smooth');

    title(' Elastic Modulus: FIT','fontsize',12)
    axis([ min_d_plot  max_d_plot -1.2 2.2])	

    text(1e-9,0.6, ['Range Indentation :' num2str(e_Indentation_IU(1)*1e9) ' to ' num2str(e_Indentation_IU(end)*1e9) '  nm'],'fontsize',12);
    text(1e-9,0.4, ['A0 :' num2str(A0*1e9) '  nm'],'fontsize',12);
    text(1e-9,0.2, ['F adhesion:' num2str(F_adhesion), ' N'],'fontsize',12)

    %%%%%% Plot Fit from Analytic 

    Force_Analytic_FIT=4/3*Chauvenete_ElasticM_mean*(R^(0.5))*(e_Indentation_IU.^(1.5));
    plot(-e_Indentation_IU,Force_Analytic_FIT/abs(F_adhesion)-1,'om','Markersize',M_size_E, 'displayname',' Analytic FIT From ElasticM (DMT): FORCE');

    text(1e-9,-0.4, ['ElasticM (DMT) AL:' num2str(Chauvenete_ElasticM_mean*1e-9) '  GPa'],'fontsize',12);
    text(1e-9,-0.6, ['ElasticM (DMT) AL (std) :' num2str(Chauvenete_ElasticM_std*1e-9) '  GPa'],'fontsize',12);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Calculate errors %%%%%%%%%%%%%%



    figure (count_figures+103)
    hold on 
    box on
    Force_Analytic_FIT_total_norm=Force_Analytic_FIT/abs(F_adhesion)-1;
    Force_total_norm =e_Force_IU(EE1:EE2)/abs(F_adhesion);
    % Difference_squared_force_fit=(((Force_total_norm-Force_Analytic_FIT_total_norm)./Force_total_norm).^2);
    Difference_squared_force_fit=(((Force_total_norm-Force_Analytic_FIT_total_norm)/max(abs(Force_total_norm))).^2);

    Data_input=Difference_squared_force_fit;

    Chauvenete;
    Sum_squared_force_fit=sum(Data_input);
    force_fit_points=length(Data_input);

    Standard_deviation_error_force_fit=floor(sqrt(1/force_fit_points*(Sum_squared_force_fit))*100);

    text(0,1.6, ['ERROR ANALTIC FIT (STD): ' num2str(Standard_deviation_error_force_fit), ' %'],'fontsize',15);


    if NonLinear_ElasticM==1

        Force_nl_reconstructed_net_norm=Force_nl_reconstructed/abs(F_adhesion)-1;
    % 	Difference_squared_force_nl_fit=(((Force_total_norm-Force_nl_reconstructed_net_norm)./Force_total_norm).^2);
        Difference_squared_force_nl_fit=(((Force_total_norm-Force_nl_reconstructed_net_norm)/max(abs(Force_total_norm))).^2);

        Data_input=Difference_squared_force_nl_fit;

        Chauvenete;
        Sum_squared_force_nl_fit=sum(Data_input);
        force_nl_fit_points=length(Data_input);

        Standard_deviation_error_force_nl_fit=floor(sqrt(1/force_nl_fit_points*(Sum_squared_force_nl_fit))*100);

        text(0,2, ['ERROR NL FIT (STD): ' num2str(Standard_deviation_error_force_nl_fit), ' %'],'fontsize',15);


    end

    %%End errors

end