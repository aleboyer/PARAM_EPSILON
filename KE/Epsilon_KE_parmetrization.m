function Epsilon_KE_parmetrization(mooring,Environment)

% mooring is a structure similar to the one find in GMACMD
% Environement

% %% set path. TODO move this when function is called
% addpath(genpath([Environment.gmacmdroot '/GMACMD/TPXO/']));
% addpath(genpath([Environment.gmacmdroot '/GMACMD/GMACMDlib/']));

%% make a time axis
mooring.timeaxis = mooring.begintime+mooring.increment/1440*(0:length(mooring.u)-1);
%% get Coriolis frequency 
speed_earth=2*pi/86400;

mooring.f=2*speed_earth*sind(-mooring.latitude); % -latitude because we are looking at rotary spectra. f is clockwise/omega negative in the northern hemisphere.
mooring.f=mooring.f*86400/2/pi;

% get TPXO tidal output to remove the barotropic tides from the recorded velocities
%   TODO: function for TPXO
% [u,~]=tmd_tide_pred([Environment.gmacmdroot,'GMACMD/TPXO/TMD/DATA/Model_tpxo7.2'],...
%     mooring.timeaxis,mooring.latitude,mooring.longitude,'u');
% [v,~]=tmd_tide_pred([Environment.gmacmdroot,'GMACMD/TPXO/TMD/DATA/Model_tpxo7.2'],...
%     mooring.timeaxis,mooring.latitude,mooring.longitude,'v');

%% stratification

% mooring  = add_n2levitus(mooring,Environment);
% mooring  = correct_timeserie(mooring); % check for plateau spikes and holes.Build blocks of more than 30 days
if isfield(mooring,'block')
    mooring  = sliding_spectra(mooring);
    mooring  = get_KE_peak(mooring,Environment);
    mooring  = get_GM(mooring,Environment);
    save(sprintf('%sGMACMD/GMACMDlocal/%s.mat',Environment.gmacmdroot,Environment.name),'mooring')
    if Environment.dsp.flag==1
        Environment.nbblock=1;
        func_GMACMDmovie(mooring,Environment);
    end
else
    Environment.noblock_ind=Environment.noblock_ind+1;
    Environment.noblock_files{Environment.noblock_ind}=Environment.file;
end


end
