Nonparametric estimation of the Fisher information

For all examples, the inputs are assumed to follow: a) mode=0, a standard Gaussian distribution; b) mode=1, a binary distribution (= 1 and -1 with probability 0.5).

[main_plot_estimates.m] is the main file to plot the functions (pdf, its derivative, the Fisher information, and the MMSE) and their estimates; 

[main_plot_complexity.m] is the main file to compare the complexities (the sample size needed to guarantee a given error and a given confidence) of the Bhattacharya estimator and those of the clipped estimator; and 

[main_err_histogram.m] is the main file to plot the bias of the two estimators.
