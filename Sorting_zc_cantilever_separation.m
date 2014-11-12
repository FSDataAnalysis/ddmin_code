%%%% Does not affect raw data

Dumb=ZsnrEx;

    if ZsnrEx(1)<0   %% First value smaller than zero 
        if ZsnrEx(1)<ZsnrEx(length(ZsnrEx))
            Dumb=-Dumb; % if increaseing zc as cantilever goes down then no more changes needed
        end
        % if decreasing zc as cantilever goes down
        
        if ZsnrEx(1)>ZsnrEx(length(ZsnrEx)) %  increaseing zc as cantilever goes down                     
           Dumb=Dumb-Dumb(end); 
        end
        
    end
    
    if ((ZsnrEx(1))>=0)   %% First value larger than zero 
        
        if ZsnrEx(1)<ZsnrEx(length(ZsnrEx)) %  increaseing zc as cantilever goes down          
            difference=2*(ZsnrEx(2:end)-ZsnrEx(1:end-1));
            for iii=2:1:length(ZsnrEx)
                Dummy=0;
                for nnn=1:1:(iii-1)
                    Dummy=Dummy-difference(nnn);
                end
                Dumb(iii)=Dumb(iii)+Dummy;
            end
            
        end
    end
      
    ZsnrEx=Dumb;