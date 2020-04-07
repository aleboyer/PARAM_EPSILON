function [E0,err]=find_Snorm(Conti,GM)    
    % compute the GM spectrum for all the E0 "guess" and compute the least square 
    Conti=Conti(:);GM=GM(:);
    GM(GM==Inf)=nan;
    inan=(~isnan(GM) & ~isnan(Conti));
    
    E0=nanmedian(GM(inan).\Conti(inan));
    err=abs(sum(Conti(inan)-E0*GM(inan)));
    
        
        
    
    

