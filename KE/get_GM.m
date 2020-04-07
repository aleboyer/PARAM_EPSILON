function mooring = get_GM(mooring,Environment)

% GM parameters
params=Gm76Params;
quant='Vel';
% BE VERY AWARE mooring.f <0 for Northern hemisphere. It was simpler to do
% this because I am looking at CW and CCW component. CW is <0 frequencies.

%% start GM fit
N             = sqrt(mooring.N2);%rad/s
f             = abs(mooring.f)*2*pi/86400;%rad/s
om=linspace(f,N,500); %rad/s
S_GM=GmOm(quant,om,f,N,params);
S_GM=S_GM(2:end)/86400*2*pi/2;om=om(2:end)*86400/2/pi;%cpd

for ii=1:length(mooring.block)
    tic
    [nbscan,~]=size(mooring.block{ii}.spec);
    omega         = mooring.block{ii}.freq;
    if mooring.f>=0 % if f<0 flip the spectra to work on a positive axis. Easier cause I do not have to devellop a negative and a positive case
        continuum = mooring.block{ii}.spec;
        edge1 = [mooring.block{ii}.KESource.edge1];
        edge2 = [mooring.block{ii}.KESource.edge2];
    else
        continuum = fliplr(mooring.block{ii}.spec);
        edge1 = fliplr(-[mooring.block{ii}.KESource.edge1]);
        edge2 = fliplr(-[mooring.block{ii}.KESource.edge2]);
    end
    
    nanomega=omega*0+1;% define a mask to nan all the source in the continuum
    for e=1:length(edge1)
        nanomega(omega>edge2(e) & omega<edge1(e))=nan;
    end
    %nanomega(abs(omega)<abs(mooring.f))=nan; % get rid of frequencies lower than f to compute the GM best fit
    
    nancontinuum = continuum; % create a continuum array that we re going to nan out
    nancontinuum(:,isnan(nanomega))=nan; % nan out the sources in the continuum
    
    nancontinuum  = nancontinuum(:,omega>0);
    continuum    = continuum(:,omega>0);
    omega        = omega(omega>0);
    
    % compute bunch of GM spectrum
    select_e0=zeros(1,nbscan);
    select_LSE=zeros(1,nbscan);
    select_inertia_slope=zeros(1,nbscan);
    select_inertia_level=zeros(1,nbscan);
    select_sub1_slope=zeros(1,nbscan);
    select_sub1_level=zeros(1,nbscan);
    select_sup1_slope=zeros(1,nbscan);
    select_sup1_level=zeros(1,nbscan);
    select_supN_slope=zeros(1,nbscan);
    select_supN_level=zeros(1,nbscan);
    % LSE is the quantity we want to minize.
    % This is the least square error between one spectrum and GM(E0)
    % eps is our futur E0
    
    for s=1:nbscan %loop over every spectra to compute GM
        Conti=nancontinuum(s,:);
        Conti(omega>10)=nan; % remove weird spectra tail wich weigh too much in the fitting
        S_GMi=interp1(om,S_GM,omega);
        Snorm_GM=S_GMi.';Snorm_GM(Snorm_GM==0)=nan;
        [E0,LSE]=find_Snorm(Conti,Snorm_GM);
        [inertia,sub1,sup1,supN] = get_slope(omega,Conti,abs(mooring.f),N/2/pi*86400);
        S1=E0*Snorm_GM;
        
        
        
        if E0==0
            fprintf('pb at %i %i\n',ii,s);
        end
        
        if Environment.dsp.flag==1
            figure(2)
            cla
            loglog(omega,continuum(s,:),'r')
            hold on
            loglog(omega,S1,'k');
            loglog(omega,10.^(inertia.poly(1)*log10(omega)+inertia.poly(2)),'m')
            loglog(omega,10.^(sub1.poly(1)*log10(omega)+sub1.poly(2)),'k')
            loglog(omega,10.^(sup1.poly(1)*log10(omega)+sup1.poly(2)),'c')
            loglog(omega,10.^(supN.poly(1)*log10(omega)+supN.poly(2)),'g')
            set(gca,'Xscale','log','Yscale','log')
            ylim([1e-8 1e0])
            pause(.05)
        end
        S1=interp1(omega,S1,mooring.block{ii}.freq);
        if mooring.f<0
            S1=fliplr(S1);
        end
        select_e0(s)=E0;select_LSE(s)=LSE;
        select_inertia_slope(s)=inertia.poly(1);select_inertia_level(s)=inertia.poly(2);
        select_sub1_slope(s)=sub1.poly(1);select_sub1_level(s)=sub1.poly(2);
        select_sup1_slope(s)=sup1.poly(1);select_sup1_level(s)=sup1.poly(2);
        select_supN_slope(s)=supN.poly(1);select_supN_level(s)=supN.poly(2);
    end
    
    mooring.block{ii}.e0=select_e0;
    mooring.block{ii}.LSE=select_LSE;
    mooring.block{ii}.inertia_slope=select_inertia_slope;
    mooring.block{ii}.inertia_level=select_inertia_level;
    mooring.block{ii}.sub1_slope=select_sub1_slope;
    mooring.block{ii}.sub1_level=select_sub1_level;
    mooring.block{ii}.sup1_slope=select_sup1_slope;
    mooring.block{ii}.sup1_level=select_sup1_level;
    mooring.block{ii}.supN_slope=select_supN_slope;
    mooring.block{ii}.supN_level=select_supN_level;
    toc
end

end

