%%%%%%%%%%%%%%%%
clc;
clear;
%%%%%%%%%%%%%%%%

% % Añadimos los path que puede necesitar el programa %%

path('Data',path)
path('Modelos',path)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Load the training model
netTrain='Model_0.mat';
load (netTrain);  

 %%Extract the dimensions of input images and the classes considered
inputSize=net.Layers(1).InputSize(1:2);
classes=net.Layers(end).Classes; 

%%Load indexes for cough images and for images with highest activation
load("CaughIndex.mat")
load ('MaxIndex.mat'); 
%%Load the spectrograms
load("SamplePatientSpectrograms.mat")

%Extract the dimmension 
num_classes = length(classes);
img_rows = inputSize(1);
img_cols = inputSize(2);

%%Load network predictions
load("Predictions_samplePatient_model0.mat")
map_patient=0;

for i=1:length(index_max)
    num_image=i;
    index=index_max(num_image);
%%Take the corresponding image
    image=spectrograms(index,:);
    image=reshape(image,img_rows,img_cols);
    image=flipud(image);

    [YPred,scores] = classify(net,image);
%%Extract the classes
    [~,topIdx] = maxk(scores, 3);
    topScores = scores(topIdx);
    topClasses = classes(topIdx);
   
    % Compute dimensions for the occlusion map
    occ_size   = 5;    % Occlusion window size (5x5 pixels)
    occ_stride = 5;    % Step size for sliding window (5 pixels per step)
    occ_pixel  = 0.5;  % Pixel value used for occlusion (0.5 for normalized image)

    map=occlusionMap(net, image,YPred, occ_size, occ_stride, occ_pixel);
    map_patient=map_patient+map;
end

% Compute the average representative map for the patient
map_patient = map_patient / length(index_max);

%Visualization of the map
figure;imagesc([0 1], [0 4.4],map_patient),colormap(jet);
set(gca,'YDir','normal'); 
ylabel('Frequency (KHz)')
xlabel('Time (s)')
colorbar;

save("MediumOcclusionMap_samplePatient", "map_patient");

