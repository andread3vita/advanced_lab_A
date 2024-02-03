### Dark Count Analysis

In this file we explain how to use the dark_count_analysis.cc script. The script was made to count and plot as a function of the applied voltage the number of dark counts in all SiPM involved in the measures.

The underlying idea it to obtain a linear dependence between the voltage and the dark counts, but the main challenge is to count properly the number of dark counts in our .root files. Writing a naive script manually is not enough because of the noise, which can be detected easily as a signal. To improve the peak identification in the histograms two functions of the TSpectrum root class were involved:

- __SmoothMarkov()__: to smooth the histograms and reduce the noise. This has pros and cons, on the one hand less noise and easier peak detection, on the other hand very close peaks can be merged together as one bigger peak.
- __Search()__: to find and count the number of peaks

The "main" function called DarkCountAnalysis looks inside the provided directory for files with name in the format: count_A_330.root where 330 is related to 33.0 V (applied voltage), and A stands for the detector (A, B, C). Than loop over all these files to obtained the desired graph of dark counts.


### Best configuration for temperature study

Since the SiPMs are quite noisy, instead of using the TSpectrum class to search for the peaks, the best strategy is to make the histograms as smooth as possible, at which point, at the cost of losing the peaks that are very close to each other, it is possible to use the algorithm manual with a threshold close to 0.

Configuration SiPM A:
- n_A = 10
- n_B = 10
- manual = true
- threshold = 0.0001
- n_search = 5

Configuration SiPM B:
- n_A = 10 
- n_B = 10
- manual = true  
- threshold = 0.005
- n_search = 5

Configuration SiPM C (both sx and dx):
- n_A = 5 
- n_B = 5
- threshold = 0.00001
- manual search
- n_search = 5
