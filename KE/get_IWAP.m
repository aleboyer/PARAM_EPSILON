%        latitude: 32.7010
%       longitude: 231.9380
%          idepth: 77
%         sfdepth: 4524
%       begintime: 727092
%         endtime: 7.2748e+05
%       increment: 60.0000
%        duration: 387.7917
%               u: [9308×1 double]
%               v: [9308×1 double]
%     temperature: [9308×1 double]
%       direction: [9308×1 double]
%        pressure: [9308×1 double]
%           speed: [9308×1 double]


Environment.gmacmdroot='/Volumes/DataDrive/';
Environment.N0=5.2e-3;
Environment.dsp.flag=0;
load([Environment.gmacmdroot 'GMACMD/Levitus/levitus_climatology.mat'],'climatology');
Environment.climatology=climatology;
load('/Volumes/DataDrive/IWAP/MP1_iwap06.mat');

Z_MP=numel(MP.z);
mooring.latitude  = nanmean(MP.lat);
mooring.longitude = nanmean(MP.lon);
mooring.idepth    = MP.z(floor(Z_MP/2));
mooring.sfdepth   = nan
mooring.begintime: 727092
mooring.endtime: 7.2748e+05
mooring.increment: 60.0000

mooring.duration: 387.7917
              mooring.u: [9308×1 double]
              mooring.v: [9308×1 double]
    mooring.temperature: [9308×1 double]
      mooring.direction: [9308×1 double]
       mooring.pressure: [9308×1 double]
          mooring.speed: [9308×1 double]



for f=1:length(Source.Source_mat)
    mooring=load([Environment.gmacmdroot Source.Source_mat{f}]);
    strname=strsplit(Source.Source_mat{f},'/');
    Environment.name=strname{end}(1:end-4);
    Environment.moviename=sprintf('%sGMACMD/MOVIE/%s.avi',Environment.gmacmdroot,Environment.name);
    Epsilon_KE_parmetrization(mooring,Environment)
    if mod(f,10)==0
        Environment.dsp.flag=1;
    else
        Environment.dsp.flag=0;
    end
end
