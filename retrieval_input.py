# details about system
#Mike's Note: 
#We recommend either fixing the stellar radius and planetary mass to the measured values, 
# or only allowing them to vary 2 standard deviations away
Target = 'HD209458'
SubFoldExtra = None #To specify strings to the end of the general subfolder name where all results are placed
Rp = 1.38  #in Jupiter radii
RpBnds = [Rp-(Rp*.2), Rp+(Rp*.2)] #+/- 20% of measured radius
Mp = 0.69 # in Jupiter mass
MpErr = None #0.03 #Error bars of mass of planet -/+, assuming symetric. From literature
T_eq =  1200.0#K
T_eqBnds = [T_eq-(T_eq*.5), T_eq+(T_eq*.5)] #+/- 50% of measured temp
R_star = 1.203  #in Solar radii
R_starErr = None #0.05 #Error bars of radius of star -/+, assuming symetric. From literature
T_star = 6092.0 #K
TsErr = None #150 #Error bars of Temp of star -/+, assuming symetric. From literature (+50K)

#parameter inital guesses and bounds. In form initia guess, then [lower bound, upper bound]
#If want to set the parameter fixed, set the lower and upper bound variable to "None"
nlive = 1000 #default is 1000
log_scatt_factor= .5
log_scatt_factorBnds = [0,1] 
logZ = 0 #log metallicity
logZBnds = [-1,3]
log_cloudtop_P = 1
log_cloudtop_PBnds = [-3,3]
CO_ratio = 0.53
CO_ratioBnds =[0.05,2.0]
scatt_slope = 1
scatt_slopeBnds = None #[-4,10]
error_multiple = 1
error_multipleBnds = [0.5, 5]
T_spot = T_star
T_spotBnds = None #[T_star-400, T_star+400] #don't expect F type star to have activity temps greater than DeltaT = 300K
spot_cov_frac = 0.0 #0.4
spot_cov_fracBnds = None #[0,0.6]

#Offsets per each insturment you have. CAN take arbitrary number of offsets. Include offset parameters even,
#even if specific insturment doesn't need offset. It needs to be included for formating purposes. 
#Just set all terms to 0 if no offsets are needed. Have offsets orded in same order of .dat file
NumInst = 3 #Should reflect number of offsets. None if absolutely no offset needed
Offset1 = 0.0001  #For IMACS/Magellan, initial offset guess
Offset1_Rng = 0.5 # fraction of MeanDepth for prior on offset
Offset1_Wavs = [3.0e-7, 1.09e-6] #Wavlength range, in meters, which insturment observes
Offset2 = 0 
Offset2_Rng = None #No offsets. Needs to be 'None' if that insturment is 'No' in .dat file
Offset2_Wavs = [1.1e-6, 1.7e-6] # for HST/WFC3, initial offset guess
Offset3 = 0 
Offset3_Rng = None #No offsets
Offset3_Wavs = [3.0e-6, 6.7e-6] #for Spitzer, initial offset guess
dataset = "HD209458b_OfSt"
