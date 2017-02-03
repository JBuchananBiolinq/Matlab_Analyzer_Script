function LOLconc = LimitOfLinearity(concentrationlist,currentlist)
    
    linear_fit = polyfit(concentrationlist(1:3),currentlist(1:3),1);
    
    MARD = zeros(length(concentrationlist));
    
    for i=1:length(concentrationlist)
        MARD(i) = 1/i*sum(abs((linear_fit(1)*concentrationlist(1:i)+linear_fit(2) - currentlist(1:i))./(linear_fit(1)*concentrationlist(1:i)+linear_fit(2))));
        if MARD(i) > 0.2 && i > 1
            x1=concentrationlist(i-1);
            y1=MARD(i-1);
            x2=concentrationlist(i);
            y2=MARD(i);
            
            LOLconc = (0.2-y1) *(x2-x1)/(y2-y1)+x1;
            break
            
        elseif i==length(concentrationlist)
            LOLconc = concentrationlist(end);
        end
    end
end
