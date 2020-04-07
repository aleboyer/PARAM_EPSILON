Environment.noblock_ind=0;
for m=1:length(Source.Source_mat)
    Environment.file=[Environment.gmacmdroot Source.Source_mat{m}];
    mooring=load(Environment.file);
%     if mooring.latitude>20 && mooring.latitude<50
    if mooring.idepth/mooring.sfdepth>.9
        break;
    end
end
