# mathematical_model_of_hydraulic_channel
Code in this repository is written in Matlab and allows to simulate passing through hydraulic channel for different types of input signals. Further it's possible to restore input signal from the output in case of noise appearence and research relationship of RMS from noise/signal ratio for different kinds of input signals.

Scripts for Matlab language:
-signal_restoring.m - main file, it contains optional filters and restores diff signals from noisy output
-RMS_err_calc.m - for calculation RMS error in case of restoring from noisy output signal
-W_line.m, bandpass_filter.m, MA_filter.m, signal_types.m - functions (see description inside of each)
