% default 0  --------------------------------------------------------------------------------------------------------%%
  
if Simulation==1

	Set_fixed_sorted=[ZsnrEx, AmEx, PhEx, DflEx, Edis_T];
     
    %%%%  Removing repeated zc values %%%%%%%%%%%%%%%%%% This is optional
    
            mmm=length(Set_fixed_sorted(:,1));
            iii=1;
            counter_removed=0;
            while iii<(mmm-1)
                ccc=iii+1;
                while ccc<(mmm-1)
                    if (Set_fixed_sorted(ccc,1))>=(Set_fixed_sorted(iii,1))
                        Set_fixed_sorted(ccc,:)=[];
                        Set_fixed_sorted(ccc,1);   %% this simply shows value on screen if double colon is removed
                        counter_removed=counter_removed+1;
                        mmm=length(Set_fixed_sorted(:,1));
                    end
                ccc=ccc+1;
                end    
                iii=iii+1;
            end
            
            ZsnrEx=Set_fixed_sorted(:,1);
            AmEx=Set_fixed_sorted(:,2);
            PhEx=Set_fixed_sorted(:,3);
            DflEx=Set_fixed_sorted(:,4);    
			Edis_T=Set_fixed_sorted(:,5); 

else


    Set_fixed_sorted=[ZsnrEx, AmEx, PhEx, DflEx];
     
    %%%%  Removing repeated zc values %%%%%%%%%%%%%%%%%% This is optional
    
            mmm=length(Set_fixed_sorted(:,1));
            iii=1;
            counter_removed=0;
            while iii<(mmm-1)
                ccc=iii+1;
                while ccc<(mmm-1)
                    if (Set_fixed_sorted(ccc,1))>=(Set_fixed_sorted(iii,1))
                        Set_fixed_sorted(ccc,:)=[];
                        Set_fixed_sorted(ccc,1);   %% this simply shows value on screen if double colon is removed
                        counter_removed=counter_removed+1;
                        mmm=length(Set_fixed_sorted(:,1));
                    end
                ccc=ccc+1;
                end    
                iii=iii+1;
            end
            
            ZsnrEx=Set_fixed_sorted(:,1);
            AmEx=Set_fixed_sorted(:,2);
            PhEx=Set_fixed_sorted(:,3);
            DflEx=Set_fixed_sorted(:,4);
end