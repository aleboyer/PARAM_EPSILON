function mooring=add_n2levitus(mooring,Environment)


% find the closest grid location to the mooring
% get the grid loc and compute the distance of each point from the mooring
% find the 4 closest grid points and make an average profile using these
% points. The average is weighted by the distance to the mooring
[LLon,LLat]=meshgrid(Environment.climatology.lon,Environment.climatology.lat);
LLon=LLon.';LLat=LLat.';
n2_avg=Environment.climatology.n2;
[LZ,LX,LY]=size(n2_avg);
n2_avg=reshape(n2_avg,[LZ,LY*LX]);
dist= sqrt((LLon-mooring.longitude).^2+(LLat-mooring.latitude).^2);
vec_dist=dist(:);
[Sort_dist,I]=sort(vec_dist);
dist_tot=sum(Sort_dist(1:4));
weight=Sort_dist(1:4)./dist_tot;
weight=ones(LZ,1)*weight.';
pr_profile=Environment.climatology.z(:,1);
mooring.N0=Environment.N0; % define N0 as the GM N0

profiles=n2_avg(:,I(1:4));
profile_n2= nansum(profiles .* weight,2);
mooring.N2=interp1(pr_profile,profile_n2,mooring.idepth);

if mooring.N2==0
    mooring.N2=nanmean(profiles(:));
    mooring.flagN2=1;
end
%% deal with surface and bottom mooring deeper and shallower than the Climato 
if(mooring.idepth>pr_profile(end)); mooring.N2=nanmin(profile_n2);end
if(mooring.idepth<pr_profile(1)); mooring.N2=profile_n2(1);end
if(isnan(mooring.N2));error('n2 is nan. Problem with interpolation');end

mooring.N2=abs(mooring.N2); % for weird point below ice at high latitude
mooring.N2bar=nanmean(profile_n2);

