% GMACMD movie
function func_GMACMDmovie(mooring,Environment)


params=Gm76Params;
quant='Vel';

%% multi taper K=3 for 95% confidence interval
ntapers=3;

close all
cmap=colormap(jet(100));
N = sqrt(mooring.N2);
f = abs(mooring.f)*2*pi/86400;%rad/s
%om=linspace(f,N,1000); %rad/s
om=logspace(log10(f),log10(N),500);
S_GM=GmOm(quant,om,f,N,params);
S_GM=S_GM(2:end)/86400*2*pi/2;om=om(2:end)*86400/2/pi;%cpd
f=mooring.f;
close all
ii=Environment.nbblock;

switch sign(f)
    case 1
        signf=1;
    case -1
        signf=-1;
    case 0
        signf=1;
end


e0=mooring.block{ii}.e0;
%% rename the variable to make it easy to write :)
freq      = mooring.block{ii}.freq;
data      = mooring.block{ii}.spec;
timeU     = datenum2yday(mooring.block{ii}.time);
time      = datenum2yday(mooring.block{ii}.timespec);

[spl,spu] = SmoothSpecCLim(freq,ntapers,0);
if f<0
    data = fliplr(data);
end

%% peaks definition
selected_loc = [mooring.block{ii}.KESource(:).fp];
edge1        = [mooring.block{ii}.KESource(:).edge1];
edge2        = [mooring.block{ii}.KESource(:).edge2];
f_loc        = find(abs([mooring.block{ii}.KESource(:).fp]-f)== ...
    min(abs([mooring.block{ii}.KESource(:).fp]-f)));
M2_loc       = find((abs(signf*[mooring.block{ii}.KESource(:).fp]-2))== ...
    min(abs(signf*[mooring.block{ii}.KESource(:).fp]-2)));

%% plot tricks
nof=0;noM2=0;
if (strcmp(mooring.block{ii}.KESource(M2_loc).name,'M2') || ...
        strcmp(mooring.block{ii}.KESource(M2_loc).name,'S2') || ...
        strcmp(mooring.block{ii}.KESource(M2_loc).name,'N2') || ...
        strcmp(mooring.block{ii}.KESource(M2_loc).name,'-M2') || ...
        strcmp(mooring.block{ii}.KESource(M2_loc).name,'-S2') || ...
        strcmp(mooring.block{ii}.KESource(M2_loc).name,'-N2'))
    b1=max(mooring.block{ii}.KESource(M2_loc).KE);
else
    b1=nan;
    noM2=1;
end
if (strcmp(mooring.block{ii}.KESource(f_loc).name,'f'))
    b2=max(mooring.block{ii}.KESource(f_loc).KE);
else
    b2=nan;
    nof=1;
end

xmax1 = max([b1 b2 ...
    max(data(:,find(freq>0,1,'first'))*1440/mooring.increment)]);
scale=10.^round(log10(xmax1./[b1 b2 ...
    max(data(:,find(freq>0,1,'first'))*1440/mooring.increment)]));
xmax  = max(data(:));
xmin  = max(min(data(data>0)),1e-10);
colorscale = linspace(0,max(abs(selected_loc)),100);

%% define seasons
ydaytime   = rem(time,365);
season     = 0*ydaytime;
season(ydaytime>=351) = 4;%winter
season(ydaytime<=81)  = 4;%winter
season(ydaytime>=81 & ydaytime<=171)  = 2;%spring
season(ydaytime>=171 & ydaytime<=261) = 1;%summer
season(ydaytime>=261 & ydaytime<=351) = 3;%autumn
seasname{1}='Summer';seasname{2}='Spring';seasname{3}='Autumn';
seasname{4}='winter';

%% plot
figure('Position',[100,100,1600,2000])
%v = VideoWriter(sprintf('../figs/Time_evolution_all_freq_rot_%i_%i_%i.avi',i,ii,iii));
v = VideoWriter(Environment.moviename);
v.FrameRate=10;
open(v)
ax(1)=subplot('Position',[.1 .47 .8 .3]);
ax(2)=subplot('Position',[.1 .1 .8 .3]);
ax(3)=subplot('Position',[.1 .85 .8 .1]);



% timeseries

hold(ax(2),'on')
if nof
    [axyy,Y1,Y2]=plotyy(ax(2),time,0*time,...
        time,e0);
    Y1.Color='w';
else
    [axyy,Y1,Y2]=plotyy(ax(2),time,mooring.block{ii}.KESource(f_loc).KE*scale(2),...
        time,e0);
    Y1.Color=cmap(find(abs(selected_loc(f_loc))<=colorscale,1,'first'),:);
end
axyy(1).YColor='k';
axyy(2).YColor='r';
axyy(1).YTick=linspace(0,xmax1,10);
if all(~isnan(e0))
    axyy(2).YTick=linspace(0,max(e0),10);
end
axyy(1).YTickLabel=sprintf('%1.2f\n',linspace(0,xmax1,10));
axyy(2).YTickLabel=sprintf('%1.2f\n',linspace(0,max(e0),10)');

Y2.Color='r';Y2.LineStyle='-.';
Y1.LineWidth=2;Y2.LineWidth=2;
if noM2==0
    lM2=plot(ax(2),time,mooring.block{ii}.KESource(M2_loc).KE*scale(1),...
        'Color',cmap(find(abs(selected_loc(M2_loc))<=colorscale,1,'first'),:),'LineWidth',2);
else
    lM2=plot(ax(2),0,0,'Color','w','LineWidth',2);
end
lM=[];itt=0;lmname=struct([]);
for s=1:length(selected_loc)
    % movie of the high peak, we plot 10*KE in the peak to see
    % something
    if (~strcmp(mooring.block{ii}.KESource(s).name,'f')  && ...
            ~strcmp(mooring.block{ii}.KESource(s).name,'M2') && ...
            signf*selected_loc(s)>2)
        itt=itt+1;
        lM(itt)=plot(ax(2),time,10*mooring.block{ii}.KESource(s).KE,...
            'Color',cmap(find(abs(selected_loc(s))<=colorscale,1,'first'),:),'LineWidth',2);
        lmname{itt}=['10 ' mooring.block{ii}.KESource(s).name];
    end
    % plot the lowest frequency (30 days^{-1}) to estimate the mesoscale
    % activity
    lM(itt+1)=plot(ax(2),time,data(:,find(freq>0,1,'first'))*1440/mooring.increment*scale(3),...
        'Color','k','LineWidth',2);
    lmname{itt+1}=sprintf('%iMeso',scale(3));
    
end

%title(ax(2),seasname{season(t)})
legend([Y1,lM2,lM,Y2],{sprintf('%if',scale(2)),sprintf('%iM2',scale(1)),lmname{:},'E_0/E_{GM}'},'location','Northeast')
set(axyy,'xlim',[time(1) time(end)])
set(axyy(1),'ylim',[0 xmax1]);
if all(~isnan(e0))
    set(axyy(2),'ylim',[0 max(e0)])
end
%set(axyy(2),'Ytick',0:.5:max(e0))
xlabel('yday');ylabel(axyy(1),'Sources (m^2 s^{-2})');
ylabel(axyy(2),'E_{GM}/E0_{GM}')

hold(ax(2),'off')
ax(2).FontSize=20;
axyy(2).FontSize=20;

% spectrum
hold(ax(1),'on')
t=1;
l1=loglog(ax(1),freq,data(t,:),'k','linewidth',2);
lgm=loglog(ax(1),om,S_GM,'m-.','linewidth',1);
limGM=plot(ax(1),[10,10],[1e-6 xmax],'k--');
plot(ax(1),freq(freq>1/15 & freq<1/5),...
    spl(freq>1/15 & freq<1/5)./(10^4),'--','Color','k','linewidth',2)
plot(ax(1),freq(freq>1/15 & freq<1/5),...
    spu(freq>1/15 & freq<1/5)./(10^4),'--','Color','k','linewidth',2)
text(nanmean(freq(freq>1/15 & freq<1/5)),mean([spl(1) spu(1)])./(10^4),'95%','Parent',ax(1))
l2=loglog(ax(1),om,S_GM.*e0(t),'r-.','linewidth',2);
loglog(ax(1),[N/2/pi*86400 N/2/pi*86400],[1e-6 xmax],'g--','linewidth',2);
hold(ax(1),'off')
legend([l1,lgm,l2,limGM],{'Obs','GM','Continuum','lim GM fit'})
xlabel(ax(1),'cpd');ylabel(ax(1),'m^2s^{-2}/cpd')

set(ax(1),'xlim',[freq(find(freq>=0,1,'first')) freq(end)])
set(ax(1),'xtick',[1./15 1./10 1./7.5 1./5 1./2.5 1 2 3 4 5 6])
set(ax(1),'xticklabel',num2str([1./15 1./10 1./7.5 1./5 1./2.5 ...
    1 2 3 4 5 6]','%1.2f'))
set(ax(1),'ylim',[.8*xmin 2*xmax ])
set(ax(1),'Xscale','log','Yscale','log')
ax(1).FontSize=20;
ax(1).XTickLabelRotation=45;
title(ax(1),sprintf('Lat=%2.1f,Lon=%3.1f,iDepth=%3.0f,sfdepth=%3.0f,N=%1.3d  (%i)',...
    mooring.latitude,mooring.longitude,mooring.idepth,mooring.sfdepth,N,Environment.mooringnumber),'fontsize',15)

% time vel
%xlabel(ax(3),'yday')
hold(ax(3),'on')
ax(3).Color=(.5+.1*season(t))*[1 1 1];
indtime=find(mooring.block{ii}.time>mooring.block{ii}.timespec(t)-15 & mooring.block{ii}.time<mooring.block{ii}.timespec(t)+15);
l3=plot(ax(3),timeU(indtime),mooring.block{ii}.uL1(indtime),...
    'k','linewidth',2);
ylabel(ax(3),'m s^{-1}')
ax(3).XLim =[time(t)-15 time(t)+15];
ax(3).YLim =[-max(abs(mooring.block{ii}.uL1)) max(abs(mooring.block{ii}.uL1))];
ax(3).FontSize=20;

delete(l1)
delete(l2)
delete(l3)

for t=1:length(time)
    % get the index to plot the time serie coresponding to the
    % plotted spectrum
    hold(ax(3),'on')
    ax(3).Color=(.5+.1*season(t))*[1 1 1];
    indtime=find(mooring.block{ii}.time>mooring.block{ii}.timespec(t)-15 & mooring.block{ii}.time<mooring.block{ii}.timespec(t)+15);
    l3=plot(ax(3),timeU(indtime),mooring.block{ii}.uL1(indtime),...
        'k','linewidth',2);
    ax(3).XLim =timeU(indtime([1 end]));
    hold(ax(3),'off')
    
    hold(ax(1),'on')
    ax(1).Color=(.5+.1*season(t))*[1 1 1];
    l1=loglog(ax(1),freq,data(t,:),'k','linewidth',2);
    for s=1:length(selected_loc)
        if (~isempty(selected_loc(s)) && sign(edge1(s))==sign(edge2(s)))
            hsource{s}=fill(signf*[edge2(s) freq(freq<=edge1(s) & freq>=edge2(s)) edge1(s)],...
                [max(min(data(:)),1e-16) mooring.block{ii}.spec(t,freq<=edge1(s) & freq>=edge2(s)) ...
                max(min(data(:)),1e-16)],cmap(find(abs(selected_loc(s))<=colorscale,1,'first'),:),'parent',ax(1));
            tsource{s}=text(signf*selected_loc(s),xmax,...
                mooring.block{ii}.KESource(s).name,...
                'backgroundcolor','w','fontsize',10,'parent',ax(1));
        end
    end
%     l2=loglog(ax(1),freq,GM(t,:),'r-.','linewidth',2);
    l2=loglog(ax(1),om,S_GM*e0(t),'r-.','linewidth',1);
%     l2=loglog(ax(1),om,GM(t,:),'r-.','linewidth',2);
%     lsl1=loglog(ax(1),freq,10.^(mooring.block{1}.inertia_slope(t)*log10(freq)+mooring.block{1}.inertia_level(t)),'g--');
    freqsub1=freq(freq>0 & freq<1);
    freqsup1=freq(freq>1 & freq<N/2/pi*86400);
    freqsupN=freq(freq>N/2/pi*86400);
    specsub1=10.^(mooring.block{1}.sub1_slope(t).*log10(freqsub1)+mooring.block{1}.sub1_level(t));
    specsup1=10.^(mooring.block{1}.sup1_slope(t).*log10(freqsup1)+mooring.block{1}.sup1_level(t));
    lsl2=loglog(ax(1),freqsub1,specsub1,'Color',[.2 .2 .2 .7],'linewidth',2);
    lsl3=loglog(ax(1),freqsup1,specsup1,'Color',[.4 .4 .4 .7],'linewidth',2);
    if ~isempty(freqsupN)
        specsupN=10.^(mooring.block{1}.sub1_slope(t).*log10(freqsupN)+mooring.block{1}.sub1_level(t));
        lsl4=loglog(ax(1),freqsupN,specsupN,'Color',[.6 .6 .6 .7],'linewidth',2);
    end

    hold(ax(1),'off')
    legend(ax(1),[l1,lgm,l2,limGM],{'Obs','GM','Continuum','lim GM fit'})

    
    hold(ax(2),'on')
    timeind=plot(ax(2),time([t t]),[0 xmax1],'k--','LineWidth',2);
    timetext=text(time(t),.85*xmax,seasname{season(t)},'backgroundcolor','y','fontsize',20,'Parent',ax(2));
    hold(ax(2),'off')

    
    
    
    frame=getframe(gcf);
    writeVideo(v,frame)
    %l1.Visible='off';l2.Visible='off';
    if t<length(time)
        delete(timeind)
        delete(timetext)
        delete(l1)
        delete(l2)
        delete(l3)
        %delete(lsl1)
        delete(lsl2)
        delete(lsl3)
        if ~isempty(freqsupN)
            delete(lsl4)
        end
        for s=1:length(selected_loc)
            delete(hsource{s});
            delete(tsource{s});
        end
    end
end
close(v)

