
if Simulation==1


%%%%% Simulation %%%%%

Minimum_error_virial=MultiplierError*min(abs(Virial_Fts_norm-Percentage_Virial));

        
Virial_element=find(abs(Virial_Fts_norm-Percentage_Virial)<=abs(5*Minimum_error_virial),1, 'last'); % this is the position of the vector element
D_min_Virial_Percentage_value=d_min(Virial_element);
Virial_Percentage_value=Virial_Fts_norm(Virial_element);

figure (count_figures+60)
hold on
plot(D_min_Virial_Percentage_value,  Virial_Percentage_value,'Vk', 'Markersize',5)

%%%%%%%% find the zero value of Virial %%%%%%%%%%%%%%
Virial_Fts_nor_cut=Virial_Fts_norm(Virial_element:end);
% Virial_element_zero=find(abs(Virial_Fts_nor_cut-Percentage_To_zero)<=abs(Minimum_error_virial_zero),1, 'last'); % this is the position of the vector element
Virial_element_zero=find(Virial_Fts_nor_cut<=0,1, 'first'); % this is the position of the vector element
D_min_Virial_Percentage_value_zero=d_min(Virial_element_zero+Virial_element-1);
Virial_Percentage_value_zero=Virial_Fts_norm(Virial_element_zero+Virial_element-1);
plot(D_min_Virial_Percentage_value_zero,  Virial_Percentage_value_zero,'Vk', 'Markersize',5);

Delta_Virial_To_Zero=D_min_Virial_Percentage_value-D_min_Virial_Percentage_value_zero;


%%%%%%%% find the one value of Virial %%%%%%%%%%%%%%
Minimum_error_virial_one=MultiplierError*min((Virial_Fts_norm-1)==0);  
Virial_element_one=find(abs(Virial_Fts_norm-1)<=abs(5*Minimum_error_virial_one),1, 'first'); % this is the position of the vector element
D_min_Virial_Percentage_value_one=d_min(Virial_element_one);
Virial_Percentage_value_one=Virial_Fts_norm(Virial_element_one);

hold on
plot(D_min_Virial_Percentage_value_one,  Virial_Percentage_value_one,'Vk', 'Markersize',5);

Delta_Virial_To_One=D_min_Virial_Percentage_value_one-D_min_Virial_Percentage_value_zero;


text(0,-1, ['Virial Error: ' num2str(Percentage_Virial*100), ' %'],'fontsize',12)
text(0,-1.4, ['DisToZero: ' num2str(Delta_Virial_To_Zero*1e9), ' nm'],'fontsize',12)
text(0,-1.8, ['DisToOne: ' num2str(Delta_Virial_To_One*1e9), ' nm'],'fontsize',12)

else
	

%%%%% experimental %%%%%

 
Minimum_error_virial=MultiplierError*min(abs(Virial_Fts_d_ss_norm_crop-Percentage_Virial));

        
Virial_element=find(abs(Virial_Fts_d_ss_norm_crop-Percentage_Virial)<=abs(5*Minimum_error_virial),1, 'last'); % this is the position of the vector element
D_min_Virial_Percentage_value=D_M_ZEROED_CROP_d(Virial_element);
Virial_Percentage_value=Virial_Fts_d_ss_norm_crop(Virial_element);

figure (count_figures+60)
hold on
plot(D_min_Virial_Percentage_value,  Virial_Percentage_value,'Vk', 'Markersize',8)

%%%%%%%% find the zero value of Virial %%%%%%%%%%%%%%
Virial_Fts_nor_cut=Virial_Fts_d_ss_norm_crop(Virial_element:end);
% Virial_element_zero=find(abs(Virial_Fts_nor_cut-Percentage_To_zero)<=abs(Minimum_error_virial_zero),1, 'last'); % this is the position of the vector element
Virial_element_zero=find(Virial_Fts_nor_cut<=0,1, 'first'); % this is the position of the vector element
D_min_Virial_Percentage_value_zero=D_M_ZEROED_CROP_d(Virial_element_zero+Virial_element-1);
Virial_Percentage_value_zero=Virial_Fts_d_ss_norm_crop(Virial_element_zero+Virial_element-1);
plot(D_min_Virial_Percentage_value_zero,  Virial_Percentage_value_zero,'Vk', 'Markersize',8);

Delta_Virial_To_Zero=D_min_Virial_Percentage_value-D_min_Virial_Percentage_value_zero;


%%%%%%%% find the one value of Virial %%%%%%%%%%%%%%
Minimum_error_virial_one=MultiplierError*min((Virial_Fts_d_ss_norm_crop-1)==0);  
Virial_element_one=find(abs(Virial_Fts_d_ss_norm_crop-(1-Percentage_Virial))<=abs(10*Minimum_error_virial_one),1, 'first'); % this is the position of the vector element
D_min_Virial_Percentage_value_one=D_M_ZEROED_CROP_d(Virial_element_one);
Virial_Percentage_value_one=Virial_Fts_d_ss_norm_crop(Virial_element_one);

hold on
plot(D_min_Virial_Percentage_value_one,  Virial_Percentage_value_one,'Vk', 'Markersize',8);

Delta_Virial_To_One=D_min_Virial_Percentage_value_one-D_min_Virial_Percentage_value_zero;

text(0,-0.4, ['Virial Error: ' num2str(Percentage_Virial*100), ' %'],'fontsize',12)
text(0,-1.2, ['DisToZero: ' num2str(Delta_Virial_To_Zero*1e9), ' nm'],'fontsize',12)
text(0,-1.8, ['DisToOne: ' num2str(Delta_Virial_To_One*1e9), ' nm'],'fontsize',12)

end


%% Experimental
% D_min_Virial_Percentage=D_M_ZEROED(Virial_element);