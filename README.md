# GAM-first-derivative-functions
Functions for visualizing first derivative of splines in BRMS

Uses the method of finite differences as suggested by Gavin Simpson to estimate the first derivative of a spline

Attempts to find value ranges with a non-zero rate of change in the predicted response

deriv_plot() annotates the regression line (posterior mean) by coloring the value ranges with a non zero rate of change
deriv_plot2() plots the actual first dervitave 
species_interact_deriv() is a temporary and brute force version of deriv_plot that allows interactions with a factor. Will eventually be deprecated and deriv_plot() will be updated to enable this more efficiently
