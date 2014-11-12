
if Simulation==1
%%%%% Simulation %%%%%

      
Attractive_element=find((abs(Fts_cons_crop_ss)-Attractive_min_force)>0,1, 'first'); % this is the position of the vector element
Damping_element=find((Damping_coefficient_crop_ss-Damping_magnitude)>0,1, 'first'); % this is the position of the vector element
D_min_Attractive_value=d_min_crop(Attractive_element);
Attractive_value=Fts_cons_crop_ss(Attractive_element);
D_min_Damping_value=d_min_crop(Damping_element);
Damping_value=Damping_coefficient_crop_ss(Damping_element);

figure (count_figures+11)
hold on
plot(D_min_Attractive_value,  Attractive_value,'Vk', 'Markersize',10)
Negativee_forces=Fts_cons_crop_ss-Error_Attractive_zero; % ;  Fts_cons_crop_ss_min
Fts_cons_crop_ss_min_element=find(Negativee_forces<Fts_cons_crop_ss_min, 1, 'first');
Adhesion_minimum=Fts_cons_crop_ss(Fts_cons_crop_ss_min_element);
D_Attractive_zero=d_min_crop(Fts_cons_crop_ss_min_element);
plot(D_Attractive_zero,  Adhesion_minimum,'Vk', 'Markersize',10)

Delta_attractive_force=D_min_Attractive_value-D_Attractive_zero;

%%% Damping %%%%%%%%%%
hold on
plot(D_min_Damping_value,  Damping_value,'Vk', 'Markersize',10)

Delta_Damping=D_min_Damping_value;


text(0,1.8*abs(F_adhesion), ['Damping ' num2str(Damping_value),' Ns/m at ' num2str(D_min_Damping_value*1e9),' nm'],'fontsize',12)
text(0,1.4*abs(F_adhesion), ['Distance attr. ' num2str(Delta_attractive_force*1e9), ' nm from ' num2str(Attractive_value*1e9), ' to ' num2str(Adhesion_minimum*1e9), ' nN'],'fontsize',12)

saveas(count_figures+11, num2str(11),'fig');

end


if Simulation==0	
  
Attractive_element=find((abs(Fts_cons_crop_sorted_ss)-Attractive_min_force)>0,1, 'first'); % this is the position of the vector element
Damping_element=find((Damping_coefficient_crop_ss-Damping_magnitude)>0,1, 'first'); % this is the position of the vector element
D_min_Attractive_value=d_min_crop_sorted_zeroed(Attractive_element);
Attractive_value=Fts_cons_crop_sorted_ss(Attractive_element);
D_min_Damping_value=d_min_Damp(Damping_element);
Damping_value=Damping_coefficient_crop_ss(Damping_element);

figure (count_figures+11)
hold on
plot(D_min_Attractive_value,  Attractive_value,'Vk', 'Markersize',10)
Negativee_forces=Fts_cons_crop_sorted_ss-Error_Attractive_zero; % ;  Fts_cons_crop_ss_min
Fts_cons_crop_ss_min_element=find(Negativee_forces<F_adhesion, 1, 'first');
Adhesion_minimum=Fts_cons_crop_sorted_ss(Fts_cons_crop_ss_min_element);
D_Attractive_zero=d_min_crop_sorted_zeroed(Fts_cons_crop_ss_min_element);
plot(D_Attractive_zero,  Adhesion_minimum,'Vk', 'Markersize',10)

Delta_attractive_force=D_min_Attractive_value-D_Attractive_zero;

%%% Damping %%%%%%%%%%
hold on
plot(D_min_Damping_value,  Damping_value,'Vk', 'Markersize',10)

Delta_Damping=D_min_Damping_value;


text(0,1.8*abs(F_adhesion), ['Damping ' num2str(Damping_value),' Ns/m at ' num2str(D_min_Damping_value*1e9),' nm'],'fontsize',12)
text(0,1.4*abs(F_adhesion), ['Distance attr. ' num2str(Delta_attractive_force*1e9), ' nm from ' num2str(Attractive_value*1e9), ' to ' num2str(Adhesion_minimum*1e9), ' nN'],'fontsize',12)

saveas(count_figures+11, num2str(11),'fig');

end


%% Experimental
% D_min_Virial_Percentage=D_M_ZEROED(Virial_element);