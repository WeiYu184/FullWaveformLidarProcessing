# FullWaveformLidarProcessing

A tool for processing the full waveform satellite lidar data. Please get started from main.m.


FUNCTION:

wave_33_prep.m  - data pre-processing, including the quantization and voltage value conversion of raw data

normalisation.m - data normalisation

bgnoise.m - removing the background noise by empirical threshold method 

transmit_processing.m - transmit waveform processing

Gaussianfilter.m - smoothing the return waveform

initial_parameter_estimation.m - estimate the initial parameters for least squares

least_squares.m - least squares estimation

QualityofFit.m - accuracy and quality evaluation

feature_extraction.m - extracting the waveform characteristics
