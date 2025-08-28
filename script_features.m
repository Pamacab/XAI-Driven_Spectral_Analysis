%%%%%%%%%%%%%%%%
clc;
clear;
%%%%%%%%%%%%%%%%

path('Data',path)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Author: Patricia Amado Caballero
% Email: patricia.amado@uva.es
% Date: 2025-08-27

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

%%Load the patients set for the calculate the spectral features
patients_set= "Weighted_Spectrograms_COPD.mat";
load(patients_set);
% Load the signal frequency parameters
fs=8820;
bands_separation=1,
if bands_separation==1
    bandLim=[10 500;500 1000;1000 1500;1500 2000;2000 fs/2];
else
    bandLim=[10 fs/2];
end

for i=1:length(spectrograms)
         spectrogram=spectrograms{i};
         t=linspace(0,1,size(spectrogram,2));
         time_lim=[t(1) t(end)];
   
        [RP_feat,SpecBand_feat,SpecCent_feat,SpecCrestFac_feat,...
        SpecEn_feat,SpecFlat_feat,SpecFlux_feat,SpecKurt_feat,...
        SpecRenyiEn_feat,SpecRolloff_feat,SpecSkew_feat]=...
        compute_spectral_features(spectrogram,fs,bandLim,t,time_lim);

        if bands_separation==1
            features_spectrogram=[RP_feat;SpecBand_feat;SpecCent_feat;SpecCrestFac_feat;...
            SpecFlat_feat;SpecFlux_feat;SpecKurt_feat;...
            SpecRenyiEn_feat;SpecRolloff_feat;SpecSkew_feat];

            features(:,:,i)=features_spectrogram;
            SpecEntropy(i,:)=SpecEn_feat;

        else
            features_spectrogram=[RP_feat,SpecBand_feat,SpecCent_feat,SpecCrestFac_feat,...
            SpecEn_feat,SpecFlat_feat,SpecFlux_feat,SpecKurt_feat,...
            SpecRenyiEn_feat,SpecRolloff_feat,SpecSkew_feat];

            features(i,:)=features_spectrogram;     

        end
end

if bands_separation==1
    save(strcat("Features_Bands_",extractAfter(extractAfter(patients_set,"_"),"_")),"features","SpecEntropy");
else
    save(strcat("Features_",extractAfter(extractAfter(patients_set,"_"),"_")),"features");
end

