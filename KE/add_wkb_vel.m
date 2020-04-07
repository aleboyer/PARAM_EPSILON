function mooring=add_wkb_vel(mooring)

% need mooring.N2 from add_n2levitus(mooring,Environment)

coef_wkb=(sqrt(mooring.N2)/mooring.N0).^-(1/2); % wkb = (N/N0)^(-1/2)

for j=1:length(mooring.block)
    mooring.block{j}.wkb=coef_wkb*(mooring.block{j}.uL1+1i*mooring.block{j}.vL1);
end
