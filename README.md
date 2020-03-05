# Average_perf_metrics
Using different averaging methods to calculate performance metrics.

This code reads in WCVI Herring 2016 model outputs to get parameters (then changes them to turn it into a longer-lived Frankenfish).

It gets the BMSY-based reference points then runs an age-structured model for 100 years. In each year it calculates whether the stock is above (1) or below (0) the LRP and USR, whether the stock is above (0) or below (1) FMSY, the catch, and whether the catch is above (1) or below (0) MSY.

The population is driven by a time series of Ft

Replicates differ from each other in terms of random noise around M, and Ft.

At the end of the period it calculates average performance metrics using three methods:
1. ICES (each and every year)
2. Average over replicates
3. Average over time and replicates

To run, just source R/01_Run_models.r

Methods of calculating performance metrics can be changed at the bottom of the file.

TODO:
- Read in a different species so that parameters don't need to be changed
- Improve the way Ft is handled
- Improve the way random errors are handled (currently just noise added to M in normal space). Add recruitment noise.
- Improve the way FMSY stats are calculated. Currently runs will crash every now an then. If this happens just source again.
- Better plots
