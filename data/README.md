# Examples of data suitable for applying this method.

## File Descriptions
### CoughIndex.mat'

This MATLAB script is the main entry point for occlusion map analysis in the context of model activations. It allow us for generating and analyzing occlusion maps from the activations of a deep learning model, using spectrograms as input data.

### 'GMM_SamplePatient.mat'

This script is responsible for performing the fundamental calculation of the occlusion map. 

### 'MaxIndex.mat'

This MATLAB script implements a 2D Gaussian Mixture Model (GMM), designed to fit and visualize complex data distributions, such as those that might be derived from an occlusion map. 

### 'Model_0.mat'

This script focuses on generating masks from occlusion maps, which allows highlighting regions of interest in the original signal. This process is fundamental for interpretability, as it helps identify which parts of the input are most relevant for a specific model decision or activation. 

### 'Predictions_samplePatient_model0.mat'

Automates the process of extracting spectral features from a dataset of spectrograms, preparing the data for subsequent analysis or modeling stages.

### 'SamplePatientSpectrograms.mat'

Function responsible for calculating various spectral features from an input spectrogram. 


### 'Weighted_Spectrogram_Pt0.7.mat, Weighted_Spectrograms_ARD.mat, Weighted_Spectrograms_CANCER.mat, Weighted_Spectrograms_COPD.mat, Weighted_Spectrograms_NonCOPD.mat'


	





	




