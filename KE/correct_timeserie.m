%function [uL1,vL1,nb_block] = correct_timeserie(mooring)
function mooring = correct_timeserie(mooring,cutperiod)
%% correct u of non physical values: 
%  exact 0 or constant values over few sample (plateau shape)
%  find out nan gap longer than 6 hours, interp everything to remove nans
%  and nan out the previsouly computed nan gap longer the 3 hours
%
%  argin: u  time seri
%         timeaxis in days
%  argout: u_L1 = correct time serie
%  Arnaud Le Boyer 03/08/2017

    %L        = length(mooring.u);
    u        = mooring.u;
    U        = mooring.utrope;
    du       = diff(u);
    v        = mooring.v;
    V        = mooring.vtrope;
    dv       = diff(v);
    
    % Need more data than 30 days. Some time series are all nans or almost
    % all nans
    if sum(isnan(mooring.v))*mooring.increment/1440<30
        
    
    if nargin==1
        cutperiod=6;% cutoff period in hours
    end
    
    u(find(du==0)+1) = nan; % deal with plateau values, not physical
    v(find(dv==0)+1) = nan; % deal with plateau values, not physical
    
    u(find(abs(du)>3*nanstd(u))+1)=nan; % remove peaks
    v(find(abs(dv)>3*nanstd(v))+1)=nan; % remove peaks
    u=filloutliers(u,'center','movmedian',10);% remove peaks
    v=filloutliers(v,'center','movmedian',10);% remove peaks
    
    %remove nans at begining and end of the time seroes for u and v.
    %I make sure u and v start and end at the same "time"
    testnan=0;
    while testnan==0
        nanu=isnan(u);nanv=isnan(v);
        ind0   =max(find(nanu==0,1,'first'),find(nanv==0,1,'first'));
        indend =min(find(nanu==0,1,'last'),find(nanv==0,1,'last'));
        u      = u(ind0:indend);
        v      = v(ind0:indend);
        U      = U(ind0:indend);
        V      = V(ind0:indend);
        
        if (~isnan(u(end)) && ~isnan(v(end)) && ...
            ~isnan(u(1))   && ~isnan(v(1)))
            testnan=1;
        end
    end
    
    % cut time axis too
    time   = mooring.timeaxis(ind0:indend);
    dt     = mooring.increment; % dt in min

    % remove the nan in the u time serie
    inddata    = find(~isnan(u));
    nan_mask   = ~isnan(u);
    mask_data  = ones(1,length(u));
    dinddata   = diff(inddata);
    ind_gap  = find(dinddata*dt/60>cutperiod);   %gap <cutperiod in hours and dt should be in minutes
    % find gap > cutperiod h - hopefully not many gaps - to define the number of blocks in our next
    % analyses
    for i=1:length(ind_gap)
        mask_data(  ...
            inddata(ind_gap(i))+1:...
            inddata(ind_gap(i))+dinddata(ind_gap(i))-1)  = nan; 
    end
    % interpolated to remove the nan in the u time serie
    uL1=interp1(time(nan_mask),u(nan_mask),time,'pchip');

    
    % remove the nan in the u time serie
    inddata    = find(~isnan(v));
    nan_mask   = ~isnan(v);
    dinddata   = diff(inddata);
    ind_gap  = find(diff(inddata)*dt/60>cutperiod);   %gap <cutperiod hours dt should be in minutes
    % find gap > cutperiod - hopefully not many gaps - likely to define the number of blocks in our next
    % analyses
     for i=1:length(ind_gap)
        mask_data(  ...
            inddata(ind_gap(i))+1:...
            inddata(ind_gap(i))+dinddata(ind_gap(i))-1)  = nan; 

    end
    vL1=interp1(time(nan_mask),v(nan_mask),time);
      
    % nan the gap > cutperiodh previsouly computed
    uL1=uL1.*mask_data;
    vL1=vL1.*mask_data;
    
    if sum(isnan(uL1))~=sum(isnan(vL1))
        warning('isnanv diff isnanu  ... weird.');
    end
    
    %% get the number of block where we have more than 30 days of continous data 
    %  and get the block 
    %  ang the timeaxis
    nonan=find(~isnan(mask_data));
    defblock=find(abs(diff(nonan))>1); % defblock should be even 
   
    defblock=[0 defblock length(nonan)]; % add first and last index to define blocks    
    uL1=uL1(:);
    vL1=vL1(:);
    count=0;
    for i=2:length(defblock)
        indedge2=nonan(defblock(i));
        indedge1=nonan(defblock(i-1)+1);
        if ((indedge2-indedge1)*dt/60>=5*30*24)
            count=count+1;
            mooring.block{count}.uL1  = uL1(indedge1:indedge2)-U(indedge1:indedge2);
            mooring.block{count}.vL1  = vL1(indedge1:indedge2)-V(indedge1:indedge2);
            mooring.block{count}.time = time(indedge1:indedge2);
        end
    end
        
 end
end