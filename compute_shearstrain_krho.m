function [P_shear,P_strain,Mmax_sh,Mmax_st,Rwtot,krho_shst,krho_st]=compute_shearstrain_krho(shearn,z_sh,strain,z_st,n2,lat,m,z_bin,Rwavg)
% compute_strain_krho.m
%
% INPUT:
%
% shearn - vector (just one profile for now) of (du/dz+i*dv/dz) normalized by a SMOOTHED version of N
% z_sh - depth vector corresponding to shearn
% strain - vector of strain
% z_st - depth vector corresponding to strain, not necc same as shear
% n2 - SMOOTHED version of N^2 that's been interpolated onto SHEAR depth vector
% lat - latitude of this profile
% m - wavenumber vector to interpolate spectra onto (so if you do this with a bunch of profiles you can compare them easily), something like 0:2*pi/300:25;
% z_bin - centers of windows over which spectra are computed, typically something like 50:150:max(z), where 50 is the base of the mixed layer
%         note that windows are half-overlapping so the 150 spacing above
%         means each window is 300 m tall
% Rwavg - shear/strain ratio of size(z_bin) if you'd like to specify it (based for example on a station average), otherwise computed using shear
%         and strain spectra for each spectral window
%
%
%  output:
%
%  P_shear - matrix of shear spectra for each depth window
%  P_strain - strain spectra
%  Mmax_sh  - cutoff wavenumber (kc) used for calculation
%  Mmax_st - cutoff wavenubmer used for strain only calculation
%  Rwtot - shear/strain ratio used, computed from spectra unless specificed in input
%  krho_shst - krho calculated from both shear and strain spectra
%  krho_st  - krho calculated from strain only, for comparison
%
%
% jen, march 09



%% Shear/strain spectra estimates of diffusivity
%   loop throgh segments of water column for spectra

delz=mean(diff(z_bin));
dz_sh=mean(diff(z_sh));
nz=length(z_bin);

K0=.05*1e-4;
N0=5.2e-3;

P_shear=NaN*ones(nz,length(m)); P_strain=P_shear;

% for ifile=1:size(shearn,2);

f=abs(2*2*pi/24/3600*sin(lat*pi/180)); 
f30=  2*2*pi/24/3600*sin(30*pi/180);

for iwin=1:nz
    
    iz=find(z_sh>=(z_bin(iwin)-delz)&z_sh<=(z_bin(iwin)+delz));
    nn=sqrt(nanmean(n2(iz)));
    Jf=f.*acosh(nn./f)/f30/acosh(N0/f30);
    
    % shear spectra
    ig=find(~isnan(shearn(iz)));
    if length(ig)>10
        [~,~,Ptot,m0] = psd_jen(detrend(shearn(iz)),mean(diff(z_sh)),'t',1);
        H=sinc(m0*mean(diff(z_sh))/2/pi).^2; % compensation for first differencing
        Ptot=Ptot./H;
        Ptot_sh=interp1(m0,Ptot,m);
        P_shear(iwin,:)=Ptot_sh(:).';
    else
        P_shear(iwin,:)=NaN; Ptot_sh=NaN*ones(size(m));
    end
    
    % strain spectra
    iz=find(z_st>=(z_bin(iwin)-delz)&z_st<=(z_bin(iwin)+delz));
    ig=find(~isnan(strain(iz))); 
    if length(ig)>10
        [~,~,Ptot,m0] = psd_jen(detrend(strain(iz)),mean(diff(z_st)),'t',1);
        H=sinc(m0*mean(diff(z_st))/2/pi).^2;
        Ptot=Ptot./H; % compensation for first differencing
        Ptot_st=interp1(m0,Ptot,m);
        P_strain(iwin,:)=Ptot_st(:).';
    else
        P_strain(iwin,:)=NaN; 
        Ptot_st=NaN*ones(size(m));
    end
    
    % find cutoff wavenumber based on SHEAR spectra
    im=find(cumsum(Ptot_sh*mean(diff(m)))<0.66);
    im=1:4;
    if length(im)>1
        Mmax_sh(iwin)=max(m(im));
    else
        Mmax_sh(iwin)=NaN;
    end
    
    % shear/strain ratio
    if ~exist('Rwavg')
        Rw=nanmean(Ptot_sh(im))/nanmean(Ptot_st(im));
        Rwtot(iwin)=Rw;
        Rwtot(Rwtot<1.01) = 1.01;
        Rw(Rw<1.01) = 1.01;
    else
        Rw=Rwavg(iwin); 
        Rwtot(iwin)=Rw;
    end
    
    hRw=3*(Rw+1)/(2*sqrt(2)*Rw*sqrt(Rw-1));
    
    Pgm=3*pi*6.3e-5*1300*3/2*m.^2./(m+3*pi/1300*nn/N0).^2; % Gm shear spectra
    
    krho_shst(iwin)=K0*sum(Ptot_sh(im)).^2./sum(Pgm(im)).^2.*hRw.*Jf;
    
    % now using strain only, use strain to compute cutoff wavenumber
    % use assumed sh/st ratio of 3 (as if we had ctd only)
    im=find(cumsum(Ptot_st*mean(diff(m)))<0.22);
    im=1:4; 
    if length(im)>1
        if max(m(im)/2/pi>.2); 
            im=find(m/2/pi<.2); 
        end
        Mmax_st(iwin)=max(m(im));
        Pgm_st=pi*6.3e-5*1300*3/2*m.^2./(m+3*pi/1300*nn/N0).^2; % GM strain spectra
        Rw=3;
        h2Rw=1/6/sqrt(2)*Rw.*(Rw+1)./sqrt(Rw-1); %for strain, different for shear
        krho_st(iwin)=K0*sum(Ptot_st(im)).^2./sum(Pgm_st(im)).^2.*h2Rw.*Jf;
    else
        Mmax_st(iwin)=NaN; krho_st(iwin)=NaN;
    end
end

