% main script to Continuum processing

% get_TPXO_barotropic_tide;
% check_extract;
% create_wkb;
% create_spectra;
% compute_peak;
% compute_GM;

Environment.gmacmdroot='/Volumes/GoogleDrive/My Drive/DATA/';
Environment.N0=5.2e-3;
Environment.dsp.flag=0;
load([Environment.gmacmdroot 'GMACMD/Levitus/levitus_climatology.mat'],'climatology');
Environment.climatology=climatology;
% % Source=load([Environment.gmacmdrot '/GMACMD/GMACMDdata/IGW_Continuum.mat']);
Source=dir('/Volumes/GoogleDrive/My Drive/DATA/GMACMD/GMACMDlocal/*.mat');
%%
% addpath(genpath('/Volumes/DataDrive/GMACMD/TPXO'));
addpath('/Volumes/GoogleDrive/My Drive/TOOLBOXES/GM/GarrettMunkMatlab/');


Environment.noblock_ind=0;
% for m=1:length(Source.Source_mat)
for m=1:length(Source)
%     try
%     Environment.file=[Environment.gmacmdroot Source.Source_mat{m}];

    Environment.file=fullfile(Source(m).folder,Source(m).name);
%     mooring=load(Environment.file);
    load(Environment.file);
%     strname=strsplit(Source.Source_mat{m},'/');Source(m).name
    Environment.name=Source(m).name(1:end-4);
    Environment.moviename=sprintf('%sGMACMD/MOVIE/%s_f%i.avi',Environment.gmacmdroot,Environment.name,m);
    Environment.mooringnumber=m;
    Epsilon_KE_parmetrization(mooring,Environment)
    fprintf('mooring %i \n',Environment.mooringnumber);
%     catch
%         fprintf('issue with mooring %i \n',Environment.mooringnumber);
%     end
    close all
end

