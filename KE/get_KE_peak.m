function mooring = get_KE_peak(mooring,Environment)

%Environment.dsp is a structure 
% dsp.flag ==1 plot and save a plot showing the different KE peaks defined
% for this data set. The peak can change between blocks and moorings

% trick for conditional cellfun
iif = @(varargin) varargin{2 * find([varargin{1:2:end}], 1, 'first')}();
    

clear TC
% define tidal component TC
TC{1}.freq=24/23.93;TC{1}.name='K1';TC{2}.freq=24/25.82;TC{2}.name='O1';
TC{3}.freq=24/24.07;TC{3}.name='P1';TC{4}.freq=24/26.87;TC{4}.name='Q1';
TC{5}.freq=24/12.42;TC{5}.name='M2';TC{6}.freq=24/12.00;TC{6}.name='S2';
TC{7}.freq=24/12.66;TC{7}.name='N2';TC{8}.freq=24/11.97;TC{8}.name='K2';
TC{9}.freq=-24/23.93;TC{9}.name='-K1';TC{10}.freq=-24/25.82;TC{10}.name='-O1';
TC{11}.freq=-24/24.07;TC{11}.name='-P1';TC{12}.freq=-24/26.87;TC{12}.name='-Q1';
TC{13}.freq=-24/12.42;TC{13}.name='-M2';TC{14}.freq=-24/12.00;TC{14}.name='-S2';
TC{15}.freq=-24/12.66;TC{15}.name='-N2';TC{16}.freq=-24/11.97;TC{16}.name='-K2';

% number of tidal component
NB_TC=length(TC);

%define tidal Harmonic
ind_TC=1:NB_TC;
it=NB_TC;
for t=ind_TC
    for tt=ind_TC(ind_TC~=t)
        it=it+1;
        TC{it}.freq=TC{t}.freq+TC{tt}.freq;
        TC{it}.name=[TC{t}.name TC{tt}.name];
    end
end

TNB_TC=length(TC);
fix_TC=TC;

clear TC;
TC=fix_TC;
% get inertial frequency
it=TNB_TC;
%02/12/2019 ALB change abs(mooring.f+TC{t}.freq) to mooring.f+TC{t}.freq
for t=ind_TC
    it=it+1;TC{it}.freq=mooring.f+TC{t}.freq;TC{it}.name=['f' TC{t}.name];
end
TC{it+1}.freq=mooring.f;TC{it+1}.name='f';TC{it+2}.freq=2*mooring.f;
TC{it+2}.name='2f';NB_TC=it+2;

for ii=1:length(mooring.block)
    freq= mooring.block{ii}.freq;
    data= mooring.block{ii}.spec;
    
    % get mean and std to determine real peaks for spectra
    peak_enhance_spec=mooring.block{ii}.meanspec; 
    % smooth the second derivative of the spectra (to define exact width of the peaks)
    Fc=6/24;% 1/4 hour filter
    [b,a]=butter(3,Fc);
    ddspec = interp1(freq(2:end-1),...
        filtfilt(b,a,diff(peak_enhance_spec,2)),freq);
    % peak height: arbitrary define as > as the max of the high freqeuncy
    % spectrum (f>5). No sources after 6 cpd
    %pks=local maxima peaks, locs= indices at wich peaks are, w=width
    %and p=prominence (height)
    peak_height=max(peak_enhance_spec(freq>6));
    [pks,locs,~,~]=findpeaks(peak_enhance_spec,freq,...
        'MinPeakProminence',.5*peak_height);

    %select main peaks (CW and CCW)
    peak_select=@(x,y) (y(pks==max(pks((abs(locs-x)==min(abs(locs-x))))))); % get the closest index from the TC freq
    
    selected_loc   = cellfun(@(x) peak_select(x.freq,locs),TC,'un',0); % associate tc frequencies to indices of peaks
    selected_peaks = cellfun(@(x) peak_select(x.freq,pks),TC,'un',0); % associate tc frequencies to peaks
    % name select peak for plot purpose
    % get the maxima peaks and other major peaks
    max_peak=max(cell2mat(selected_peaks));
    % center the freq vecteur on the peak
    center_freq=@(x) (freq-x);
    fmfp=cellfun(center_freq,selected_loc,'un',0);
    ind_fp=cellfun(@(x) find(x==0),fmfp,'un',0);
    indp=cellfun(@(x) (x>0),fmfp,'un',0);
    indm=cellfun(@(x) (x<0),fmfp,'un',0);
    % take both side of the first derivative of the spectra
    % to get the edges of the peak (~ 0 of the 1st derivative)
    % WATCH OUT: we suppose that a peak is alway surrounded be smaller
    % peak and thus the change of sign in the derivative can help us to
    % define its shape (width)
    % frequency of a peak must fit into edge 1 and edge 2
    clear ddspecp; clear ddspecm; clear ind_edge1;
    clear ind_edge2; clear edge1; clear edge2;
    ddspecp=cellfun(@(x) ddspec(x),indp,'un',0);
    ddspecm=cellfun(@(x) ddspec(x),indm,'un',0);
    ind_edge1 = cellfun(@(x) find(x>0,1,'first'),ddspecp,'un',0); % always 1 but with this I am sure to start at the right place
    edge1     = cellfun(@(x,y) freq(x+y+1),ind_fp,ind_edge1,'un',0);
    ind_edge2 = cellfun(@(x) find(flipud(x')>0,1,'first'),ddspecm,'un',0);% always 1 but with this I am sure to start at the right place
    edge2     = cellfun(@(x,y) freq(x-y-1),ind_fp,ind_edge2,'un',0);
    %% if a peak is to close to the lower freq edge of the spectrum
    % when two peaks correpond to same freq the bigger is taken into
    % account only
    edge2=cellfun(@(x,y,z) iif(isempty(x),max(freq(1),2*freq(y)-z),~isempty(x),x),...
        edge2,ind_fp,edge1,'un',0);
    
    % fun with flag :) used to get the name of the peak
    create_flag=@(x,y,z) iif((x.freq<y & x.freq>z & x.freq~=0),1,...
        (x.freq>y | x.freq<z | x.freq==0),0);
    
    % if any of the theoritcal tidal componant are included in the
    % peak definition flag=1, if not flag=0. If inertia peak == one
    % of the other peak flag =0
    clear flag;
%     flag=cellfun(create_flag,TC,edge1,edge2,'un',0);
    for i=1:length(edge2)
        if (TC{i}.freq<edge1{i} & TC{i}.freq>edge2{i} & TC{i}.freq~=0)
            flag{i}=1;
        end
        if (TC{i}.freq>edge1{i} | TC{i}.freq<edge2{i} | TC{i}.freq==0)
            flag{i}=0;
        end
    end
    
    selected_loc=cell2mat(selected_loc);edge1=cell2mat(edge1);
    selected_peaks=cell2mat(selected_peaks);
    edge2=cell2mat(edge2);flag=cell2mat(flag);
    selected_loc=selected_loc(flag==1);edge1=edge1(flag==1);
    edge2=edge2(flag==1);selected_peaks=selected_peaks(flag==1);
    clear KESource;
    KESource=[TC{flag==1}];
        
    [selected_loc,IA,~]=unique(selected_loc);
    selected_peaks=selected_peaks(IA);
    edge1=edge1(IA);edge2=edge2(IA);
    KESource=KESource(IA);
    [edge1,IA,~]=unique(edge1);
    selected_peaks=selected_peaks(IA);
    selected_loc=selected_loc(IA);
    edge2=edge2(IA);KESource=KESource(IA);
    % deal with adjacent peaks. avoid count 2 times the KE of the
    % adjacent frequency
    if (any(edge1(1:end-1)-edge2(2:end))>=0 && length(edge1)>2) % 2nd condition needed because wierd behaviour of any(x) when x is a scalar
        ind_adj1=edge1(1:end-1)-edge2(2:end)>=0;
        ind_adj2=arrayfun(@(x,y) find(freq==min([x,y])),edge1(ind_adj1),edge2([false ind_adj1]));
        edge1(ind_adj1)=freq(ind_adj2-1);
    end
    % force the closest peak to f to be the f peak
    if  any(edge1>=mooring.f & edge2<=mooring.f)
        bool1=edge1>mooring.f & edge2<mooring.f;
        if (sign(edge1(bool1)) == sign(edge2(bool1)))
            KESource(edge1>mooring.f & edge2<mooring.f).name='f';
            lowf_continuum=edge2(find(edge2<mooring.f,1,'last'));
        else
            lowf_continuum=edge2(selected_peaks==max(selected_peaks));
        end
    else
        lowf_continuum=edge2(selected_peaks==max(selected_peaks));
    end
    
    % force the closest peak to 2f to be the 2f peak
    if  any(edge1>=2.*mooring.f & edge2<=2.*mooring.f)
        bool1=edge1>2.*mooring.f & edge2<2.*mooring.f;
        if (sign(edge1(bool1)) == sign(edge2(bool1)))
            KESource(edge1>2.*mooring.f & edge2<2.*mooring.f).name='2f';
        end
    end

    
    
    if Environment.dsp.flag==1
        % define colormap for plots and earth rotation speed
        figure
        cmap=colormap(jet(100));
        
        plot(freq,peak_enhance_spec,'k','linewidth',2);
        hold on
        for t=1:NB_TC
            plot([TC{t}.freq,TC{t}.freq],[min(peak_enhance_spec),max(peak_enhance_spec)],'k--','linewidth',.5)
        end
        fill([-lowf_continuum freq(freq>=-lowf_continuum & freq<=lowf_continuum) lowf_continuum],...
            [min(peak_enhance_spec) peak_enhance_spec(freq>=-lowf_continuum & freq<=lowf_continuum) ...
            min(peak_enhance_spec)],[1 1 1])
        fill([freq(1) freq(freq<=-lowf_continuum) -lowf_continuum],...
            [min(peak_enhance_spec) peak_enhance_spec(freq<=-lowf_continuum) ...
            min(peak_enhance_spec)],.2*[1 1 1])
        fill([lowf_continuum freq(freq>=lowf_continuum) freq(end)],...
            [min(peak_enhance_spec) peak_enhance_spec(freq>=lowf_continuum) ...
            min(peak_enhance_spec)],.2*[1 1 1])
        colorscale=linspace(0,max(abs(selected_loc)),100);
        for s=1:length(selected_loc)
            fill([edge2(s) freq(freq<=edge1(s) & freq>=edge2(s)) edge1(s)],...
                [min(peak_enhance_spec) peak_enhance_spec(freq<=edge1(s) & freq>=edge2(s)) ...
                min(peak_enhance_spec)],cmap(find(abs(selected_loc(s))<=colorscale,1,'first'),:))
            text(selected_loc(s),selected_peaks(s), KESource(s).name,'backgroundcolor','y','fontsize',25)
        end
        plot(mooring.f.*[1 1],[min(peak_enhance_spec),max(peak_enhance_spec)],'k--','linewidth',2)
        hold off
        grid on
        set(gca,'xlim',[-8,8])
        set(gca,'fontsize',20)
        set(gca,'ylim',[0,1.1*max_peak])
        set(gca,'yscale','log','xscale','linear')
        xlabel('cpd','fontsize',20)
        ylabel('m^2 s^{-2} /cpd','fontsize',20)
        fig=gcf;fig.PaperPosition=[0 0 50 15];
        title(sprintf('%s,%i-block%i',Environment.name,ii))
        print([Environment.gmacmdroot 'GMACMD/FIGURE/Get_peak' Environment.name '_block' num2str(ii) '.png'],'-dpng2')
    end
    % compute the KE as function of "time" one point every 24 hours
    % (check sliding spectra)
    freq_res=nanmean(diff(freq));
    for s=1:length(selected_loc)
        if size(data,1)==1 % case where there the total time serie is 30 day (~ one spectrum)
            KESource(s).KE=sum(data(freq<=edge1(s) & freq>=edge2(s))).* freq_res;
        else % general case
            KESource(s).KE=sum(data(:,freq<=edge1(s) & freq>=edge2(s)),2).* freq_res;
        end
        KESource(s).fp=selected_loc(s);
        KESource(s).edge1=edge1(s);
        KESource(s).edge2=edge2(s);        
    end
    KESource(s+1).name='Mesoscale-CW';
    KESource(s+1).KE=sum(data(:,freq<=.066 & freq>=0.033),2).* freq_res;
    KESource(s+1).fp=0.045;
    KESource(s+1).edge1=.066;
    KESource(s+1).edge2=0.033;

    KESource(s+2).name='Mesoscale-CCW';
    KESource(s+2).KE=sum(data(:,freq>=-.066 & freq<=-0.033),2).* freq_res;
    KESource(s+2).fp=-0.045;
    KESource(s+2).edge2=-.066;
    KESource(s+2).edge1=-0.033;

    
    mooring.block{ii}.KESource=KESource;
end
       

