# mathematical_model_of_hydraulic_channel
Code in this repository is written in Matlab and allows to simulate passing through hydraulic channel for different types of input signals. Further it's possible to restore input signal from the output in case of noise appearence and research relationship of RMS from noise/signal ratio for different kinds of input signals.

Scripts for Matlab:
-input_output_redacture.m - main file, it contains optional filters and restores diff signals from noisy output
-MA_filter_testing - just attempts of using this type of filter
-W_line_characteristics.m - amplitude-freequency and phase-freequency char. for given hydraulic channel
-RMS_err_calc.m - for calculation RMS error in case of restoring from noisy output signal
