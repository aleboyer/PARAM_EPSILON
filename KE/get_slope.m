function [inertia,sub1,sup1,supN] = get_slope(omega,nancontinuum,f,N)

om_sub1=omega(omega<=1);
sub1_b=nancontinuum(omega<=1);
om_sup1=omega(omega>=3 & omega<=10);
sup1_b=nancontinuum(omega>=3 & omega<=10);


om_inertia=omega(omega>=f & omega<=N);
inertia_b=nancontinuum(omega>=f & omega<=N);



om_supN=omega(omega>=N);
if numel(om_supN)>5
    supN_b=nancontinuum(omega>=N);
    supN.poly=polyfit(log10(om_supN(~isnan(supN_b))),log10(supN_b(~isnan(supN_b))),1); 
else
    supN.poly=[nan nan];
end
inertia.poly=polyfit(log10(om_inertia(~isnan(inertia_b))),log10(inertia_b(~isnan(inertia_b))),1); 
sub1.poly=polyfit(log10(om_sub1(~isnan(sub1_b))),log10(sub1_b(~isnan(sub1_b))),1); 
sup1.poly=polyfit(log10(om_sup1(~isnan(sup1_b))),log10(sup1_b(~isnan(sup1_b))),1); 

