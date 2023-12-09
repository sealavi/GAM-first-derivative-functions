# GAM-first-derivative-functions
Functions for visualizing first derivative of splines in BRMS

Uses the method of finite differences as suggested by Gavin Simpson to estimate the first derivative of a spline


Attempts to find value ranges with a non-zero rate of change in the predicted response

deriv_plot() annotates the regression line (posterior mean) by coloring the value ranges with a non zero rate of change

deriv_plot2() plots the actual first dervitave 

species_interact_deriv() is a temporary and brute force version of deriv_plot that allows interactions with a factor. Will eventually be deprecated and deriv_plot() will be updated to enable this more efficiently

deriv_plot_zprob() Allows for estimation of derivative at different confidence levels set by the user. Z-transforms inputs to make epsillons of each dimension comparable. Provides information on how much of the posterior distrubution is on either side of zero.

deriv_ranges() extracts ranges of derivative at one side of 0

deriv_plot_by() plots the first dervitave of a one-dimensional spline with a by-variable interaction
