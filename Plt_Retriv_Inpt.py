############# General details about system #############
#Note: 
#We recommend either fixing the stellar radius and planetary mass to the measured values, 
# or only allowing them to vary 2 standard deviations aw
Rp = 1.38   #in Jupiter radii
RpBnds = [Rp-(Rp*.2), Rp+(Rp*.2)] #+/- 20% of measured radius
Mp =  0.69 # in Jupiter mass
MpErr = None #0.03 #Error bars of mass os planet -/+, assuming symetric. From literature
T_eq =  1200.0 #K
T_eqBnds = [T_eq-(T_eq*.5), T_eq+(T_eq*.5)] #+/- 50% of measured temp
R_star = 1.203  #in Solar radii
R_starErr = None #0.05 #Error bars of radius of star -/+, assuming symetric. From literature
T_star = 6092.0 #K
TsErr = None #150

############## Retrieval parameter inital guesses and bounds #############
#In form initia guess, then [lower bound, upper bound]
#If want to set the parameter fixed, set the lower and upper bound variable to "None"
nlive = 1000 #Number of live samples for Nested Sampeling. default is 1000
log_scatt_factor= .5
log_scatt_factorBnds = [0,1] 
logZ = 0 #log metallicity
logZBnds = [-1,3]
log_cloudtop_P = 1
log_cloudtop_PBnds = [-3,3]
CO_ratio = 0.53
CO_ratioBnds = [0.05,2.0]
scatt_slope = 1
scatt_slopeBnds = None #[-4,10]
error_multiple = 1 #2
error_multipleBnds = [0.5, 5]
T_spot = T_star
T_spotBnds = None #[T_star-400, T_star+400] #don't expect F type star to have activity temps greater than DeltaT = 300K
spot_cov_frac = 0.0 #0.4
spot_cov_fracBnds = None #[0,0.6]
#Offsets
InstNum = 3 #None if no offsets
Offset1 = 0.0001
Offset1_Rng = 0.5 # fraction of MeanDepth for prior on offset
Offset1_Wavs = [3.0e-7, 1.09e-6] 
Offset2 = 0 
Offset2_Rng = None
Offset2_Wavs = [1.1e-6, 1.7e-6] 
Offset3 = 0 
Offset3_Rng = None 
Offset3_Wavs = [3.0e-6, 6.7e-6]

############## Plotting setting #############
N = 1000 # number of elements for each chain to resample. For plotting only. Advised same as nlive
PlotDepth = False # to plot in terms of depth (ppm) or RpRs
X_bounds = None #Set the x and y bounds of plot. If set to None, bounds will be determined by data
Y_bounds = [0.118,0.1271] # in same units as plot units

############## Path info #############
Target = 'HD209458'
dataset = "HD209458b_OfSt"
SubFolder = 'Test'