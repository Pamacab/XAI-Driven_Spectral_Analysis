clc;
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Author: Patricia Amado Caballero
% Email: patricia.amado@uva.es
% Date: 2025-08-27

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path('Data',path)

% Load the average maps and the spectrograms for each patient
load('MediumOcclusionMap_samplePatient.mat')
load("SamplePatientSpectrograms.mat")

% Select the percentage of the signal of interest to highlight important regions, then search by patient
maskPercentage=0.7;
spectrogram=mediumSpectrogram;    
% Compute the mean signal power of the occlusion map to obtain mask thresholds
pow=bandpower(map_patient);
maxPower=max(pow);
indexes=find(pow>=maskPercentage*maxPower);
 % Identify points of maximum power in the occlusion map
PowerValue=map_patient(indexes);
minThreshold=max(PowerValue);
minIndex=find(minThreshold==map_patient);


% Compute the maximum of the spectrogram and occlusion map, then normalize
maxPatient=max(max(map_patient));
maxSpectrogram=max(max(spectrogram));
map_patientNorm=map_patient/maxPatient;
map_patientOri=map_patient;
map_patient=map_patientNorm;
maxThreshold=max(max(map_patient));
% Get the quartile containing the selected percentage of signal in the map
 aux=quantile(map_patient,maskPercentage);
 threshold2= quantile(aux,maskPercentage);
% Calculate the mask and the resulting masked spectrum
mask=map_patient>threshold2;
map_patient_m=mask.*map_patient;
spectrogramNorm=spectrogram/maxSpectrogram;
spectrogram=spectrogramNorm;
weighted_spectrogram=spectrogram.*(mask);
%Visualization
figure(1);imagesc(spectrogram);colormap(jet);
figure(2); imagesc(mask)          
figure(3); imagesc(map_patient_m)
colormap(jet)
min_val=min(min(spectrogram));
max_val=max(max(spectrogram));
figure(4); imagesc(weighted_spectrogram) ,colormap(jet); colorbar;caxis([min_val max_val])
    

name=strcat('Weighted_Spectrogram_Pt',num2str(maskPercentage),'.mat');
save(name,'weighted_spectrogram','maskPercentage');

