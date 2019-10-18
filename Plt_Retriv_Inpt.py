# details about system
#Note: 
#We recommend either fixing the stellar radius and planetary mass to the measured values, 
# or only allowing them to vary 2 standard deviations aw
Target = 'WASP31'
SubFoldExtra = None #To specify extra strings to general subfolder naming scheme
Rp = 1.549  #in Jupiter radii
RpBnds = [Rp-(Rp*.2), Rp+(Rp*.2)] #+/- 20% of measured radius
Mp = 0.478 # in Jupiter mass
MpErr = None #0.03 #Error bars of mass os planet -/+, assuming symetric. From literature
T_eq = 1575 #K
T_eqBnds = [T_eq-(T_eq*.5), T_eq+(T_eq*.5)] #+/- 50% of measured temp
R_star = 1.252  #in Solar radii
R_starErr = None #0.05 #Error bars of radius of star -/+, assuming symetric. From literature
T_star =  6300 #K
TsErr = None #150

#parameter inital guesses and bounds. In form initia guess, then [lower bound, upper bound]
#If want to set the parameter fixed, set the lower and upper bound variable to "None"
nlive = 1000 #Number of live samples for Nested Sampeling. default is 1000
N = 1000 # number of elements for each chain to resample. For plotting only. Advised same as nlive
PlotDepth = False
log_scatt_factor= .5
log_scatt_factorBnds = [-10,10] 
logZ = 0 #log metallicity
logZBnds = None #[0,3]
log_cloudtop_P = 1
log_cloudtop_PBnds = [-3,1]
CO_ratio = 0.53
CO_ratioBnds = None #[0.05,2.0]
scatt_slope = 1
scatt_slopeBnds = [-4,10]
error_multiple = 1 #2
error_multipleBnds = [0.5, 5]
T_spot = T_star
T_spotBnds = None #[T_star-400, T_star+400] #don't expect F type star to have activity temps greater than DeltaT = 300K
spot_cov_frac = 0.0 #0.4
spot_cov_fracBnds = None #[0,0.6]
#Offsets
InstNum = 3 #if no offsets
Offset1 = 0.0001 
Offset1_Rng = 0.5 # fraction of MeanDepth for prior on offset
Offset1_Wavs = [3.5e-7, 1e-6]
Offset2 = 0 
Offset2_Rng = None
Offset2_Wavs = [1.1e-6, 1.7e-6] 
Offset3 = 0 
Offset3_Rng = None 
Offset3_Wavs = [3.0e-6, 4.0e-6]
dataset = "ExoReCmbind_LC114QdEd"
