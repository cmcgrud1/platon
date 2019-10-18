import pickle
import matplotlib
from matplotlib import rc,rcParams
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import matplotlib.pyplot as plt
from dynesty import utils as dyutil
import Plt_Retriv_Inpt as PRI
import corner
from platon.transit_depth_calculator import TransitDepthCalculator
from platon.constants import R_sun, R_jup, M_jup
import seaborn as sns
import os

N = PRI.N # number of elements for each chain to resample, advised 1000
if N < 10:
	N = 10 # can't have less than 10 samples
PlotDepth = PRI.PlotDepth #if you want to plot in terms of Depth (ppm) or RpRs

#To load info that the retriever code used
direc_name = PRI.Target+'/'

#List of all compute_depths input parameters, in the right order
ParamNames = ['Rs', 'Mp', 'Rp', 'T', 'logZ', 'CO_ratio', 'log_cloudtop_P', 'log_scatt_factor', 'scatt_slope', 'T_star', 'T_spot', 'spot_cov_frac'] # error_multiple can be fit for in ln_like, but that's after computing depth 
InstNum = PRI.InstNum
Params0 =[PRI.R_star*R_sun, PRI.Mp*M_jup, PRI.Rp*R_jup, PRI.T_eq, PRI.logZ, PRI.CO_ratio, PRI.log_cloudtop_P, PRI.log_scatt_factor, PRI.scatt_slope, PRI.T_star, PRI.T_spot, PRI.spot_cov_frac]

#To add offset terms if they exist
if InstNum:
	ParamNames.append('InstNum')
	Params0.append(InstNum)
	for i in range(InstNum):
		exec("ParamNames.append('Offset"+str(i+1)+"')")
		exec("Params0.append(PRI.Offset"+str(i+1)+")")
		exec("ParamNames.append('Offset"+str(i+1)+"Dwn_Wavs')") # Need to seperate upper and lower wavelengths of insturment because need to easily put values in array
		exec("Params0.append(PRI.Offset"+str(i+1)+"_Wavs[0])")
		exec("ParamNames.append('Offset"+str(i+1)+"Up_Wavs')") #Not same name as in compute_depths, but doesn't matter because will always be fixed anyways
		exec("Params0.append(PRI.Offset"+str(i+1)+"_Wavs[1])")

#To load the data from .dat file
Wav_low, Wav_up, Depth, Err_up, Err_dwn, Insturment, Offset = np.loadtxt('Dats/'+PRI.dataset+'.dat', dtype = np.str, unpack = True)
Insturments = {}
instCnt = 1
stp = 0
MaxWav =-1
MinWav = 1e19
MaxDepth = -1
MinDepth = 1e19
bins_up =  np.array([])
bins_dwn = np.array([])
for ins in Insturment:
	#To find the bounds of the plots
	depth = float(Depth[stp])
	wav_lw, wav_up =float(Wav_low[stp]), float(Wav_up[stp])
	err_dwn, err_up= float(Err_dwn[stp]), float(Err_up[stp])
	if wav_lw < MinWav: # Used to find the bounds for plotting purposes
		MinWav = wav_lw
	if wav_up > MaxWav:
		MaxWav = wav_up
	if depth-err_dwn < MinDepth:
		MinDepth = depth-err_dwn
	if depth+err_up > MaxDepth:
		MaxDepth = depth+err_up 
	#actually storing .dat data
	bins_up = np.append(bins_up, wav_up*1e-4) # to convert from angstrom to meters
	bins_dwn = np.append(bins_dwn, wav_lw*1e-4) 
	if ins not in Insturments.keys():
		Insturments[ins] = {}
		Insturments[ins]['count'] = instCnt
		instCnt += 1
		Insturments[ins]['wav low'] = np.array(wav_lw)
		Insturments[ins]['wav up'] = np.array(wav_up)
		Insturments[ins]['Depth'] = np.array(depth)
		Insturments[ins]['Err up'] = np.array(err_up)
		Insturments[ins]['Err dwn'] = np.array(err_dwn)
	else:
		Insturments[ins]['wav low'] = np.append(Insturments[ins]['wav low'], wav_lw)
		Insturments[ins]['wav up'] = np.append(Insturments[ins]['wav up'], wav_up)
		Insturments[ins]['Depth'] = np.append(Insturments[ins]['Depth'], depth)
		Insturments[ins]['Err up'] = np.append(Insturments[ins]['Err up'], err_up)
		Insturments[ins]['Err dwn'] = np.append(Insturments[ins]['Err dwn'], err_dwn)
	stp += 1
MinWav, MaxWav = MinWav*1e-4, MaxWav*1e-4 # convert from angst to microns
bins = np.zeros((len(bins_up),2))
bins[:,0] = bins_dwn
bins[:,1] = bins_up

#List of fixed parameter. In the same order as 'ParamNames'
#If fixed, value = True. When Fixed then the corresponding Params0 value will be used for every iteration of compute_depth
#If fitted for, value = False. When not fixed the corresponding Parameter value will be drawn from the sampled chains 
ParamFixed = [True]*len(Params0)
# To get the proper directory name based on you input options
direc_name += "FIT_"
if PRI.R_starErr:
	direc_name += "Rs"
	ParamFixed[0]=False
if PRI.MpErr:
	direc_name += "Mp"
	ParamFixed[1]=False
if PRI.RpBnds:
	direc_name += "Rp"
	ParamFixed[2]=False
if PRI.T_eqBnds:
	direc_name += "Tp"
	ParamFixed[3]=False
if PRI.logZBnds:
	direc_name += "Lz"
	ParamFixed[4]=False
if PRI.CO_ratioBnds:
	direc_name += "Co"
	ParamFixed[5]=False
if PRI.log_cloudtop_PBnds:
	direc_name += "Lp"
	ParamFixed[6]=False
if PRI.log_scatt_factorBnds:
	direc_name += "Ls"
	ParamFixed[7]=False
if PRI.scatt_slopeBnds:
	direc_name += "Ss"
	ParamFixed[8]=False
if PRI.error_multipleBnds:
	direc_name += "Em"		#NOT ENTIRELY SURE HOW TO APPLY erro_multiple best fit to final plot.  
# 	ParamFixed[9]=False     #TODO: figure out how to implement error_multiple fit results to final plot
if PRI.TsErr:
	direc_name += "Ts"
	ParamFixed[9] = False
if PRI.T_spotBnds:
	direc_name += "Th"
	ParamFixed[10]=False
if PRI.spot_cov_fracBnds:
	direc_name += "Cf"
	ParamFixed[11]=False
if InstNum:
	InsCnt = 12 # To keep track of what param number 'InstNum' is, so I can add from there to the offset terms
	direc_name += 'Os'
	CNT = 1
	for n in range(InstNum):
		exec("Rng = PRI.Offset%s_Rng" %(str(n+1)))
		if Rng:
			ParamFixed[InsCnt+CNT]=False #skip InsCnt because the number of insturments is always fixed
		CNT +=3 # Do every 3rd parameter because the wavelength range (upper and lower bound) of insturment is always fixed
if PRI.nlive != 1000: # the default nlives is 1000
    direc_name += 'N'+str(PRI.nlive)
if PRI.SubFoldExtra: #To add extra string values to general subfolder naming scheme
    direc_name += PRI.SubFoldExtra+'/'
else:
    direc_name += '/'

#Loading multinest results
results = pickle.load(open(direc_name+'multinest_result.pkl', 'r'))
fit_info = pickle.load(open(direc_name+'multinest_FitInfo.pkl', 'r'))
param_names = fit_info.fit_param_names

########################To resample depths and other parameters########################
if PlotDepth: #if you want to plot in terms of Depth (ppm) or RpRs
	SaveAs = '_Depth'
else:
	SaveAs = '_RpRs'
#To check to if all plotting files have been created already. If not do resampling to create them
Results1Exist = os.path.isfile(direc_name+'results_full_outputN'+str(N)+'.pkl') #This is the number of samples for resampling DIFFERENT FROM NLIVE!
Results2Exist = os.path.isfile(direc_name+'RetrievalFits'+SaveAs+'_N'+str(N)+'.dat')
Results3Exist = os.path.isfile(direc_name+'FullResultsN'+str(N)+'.pkl')
if Results1Exist and Results2Exist and Results3Exist:
	print "Found 'results_full_outputN"+str(N)+".pkl', 'RetrievalFits"+SaveAs+"_N"+str(N)+".dat', and FullResultsN"+str(N)+".pkl file."
	print "Using '"+direc_name+"retrievalFits"+SaveAs+"_N"+str(N)+".dat' file to make plots... \n"
	wavLen,median, lower, upper = np.loadtxt(direc_name+'RetrievalFits'+SaveAs+'_N'+str(N)+'.dat', unpack = True)
	full_output = pickle.load(open(direc_name+'results_full_outputN'+str(N)+'.pkl', 'r'))
	FullResults = pickle.load(open(direc_name+'FullResultsN'+str(N)+'.pkl', 'r'))

#Only do resampling if hasn't been done for this set up
else:	
	print "Didn't find 'results_full_outputN"+str(N)+".pkl', 'RetrievalFits"+SaveAs+"_N"+str(N)+".dat', and FullResultsN"+str(N)+".pkl file."
	print "Resampling transit depth parameters based on '"+direc_name+"multinest_result.pkl' file... \n"
	#To pull out the sampled chains for each fitted for parameter
	chain = dyutil.resample_equal(results['samples'], results['weights'])

	#To get parameters in right arrangement needed for compute_depths(). i.e put all the fixed params in
	#For the values that were not fitted for set equal to the initial values (Params0)
	FullResults = np.ones((results['samples'].shape[0],len(ParamNames)))
	print "results['samples'].shape", results['samples'].shape
	print ' param_names:',  param_names
	tstCont = 0
	for p in range(len(ParamFixed)):
		if ParamFixed[p]:
			FullResults[:,p] = FullResults[:,p]*Params0[p]
		else:
			p_idx = param_names.index(ParamNames[p])
			FullResults[:,p] = chain[:,p_idx]
		tstCont +=1

	#To get the best fit model and 'confidence bars'
	calculator = TransitDepthCalculator()
	all_depths = [0]*N
	step = 0
	prntStp = N/10 #to scale the printing to always give about 10 updates
	for params in FullResults[0:N]:
		if (step%prntStp) == 0: # to print update
			print 'step:'+str(step)+'/'+str(N-1)
		cld_tp_P = 10**params[6] #covert log_cloudtop_press to cloudtop_press
		sct_fct = 10**params[7] #covert log_scatt_fract to scatt_fract
		CmpDeptTxt = "wavLen, depths, info_dict = calculator.compute_depths(params[0],params[1],params[2],params[3],logZ=params[4],"
		CmpDeptTxt += "CO_ratio=params[5],cloudtop_pressure=cld_tp_P,scattering_factor=sct_fct,scattering_slope=params[8]," # error_multiple = params[9],"
		CmpDeptTxt += "T_star=params[9], T_spot=params[10], spot_cov_frac=params[11], full_output=True"
		if InstNum: #to add offsets if need be
			CmpDeptTxt += ", InstNum = InstNum"
			cnt =1
			for I_N in range(InstNum):
				CmpDeptTxt += ', Offset'+str(I_N+1)+'=params['+str(InsCnt+cnt)+']'
				CmpDeptTxt += ', Offset'+str(I_N+1)+'_Wavs=[params['+str(InsCnt+cnt+1)+'],params['+str(InsCnt+cnt+2)+']]' 
				cnt += 3
		exec(CmpDeptTxt+')') # To compute the transit spectrum and full_output parameters for each chain
		if step == 0: #To initialize full_output dictionary
			full_output = {}
			for k in info_dict.keys():
				full_output[k] = [0]*N
		for key in full_output.keys(): #To fill up full_output dictionary
			full_output[key][step] = info_dict[key]
		all_depths[step] = depths #To save all the depth values of the resampling
		step += 1
	all_depths = np.array(all_depths)

	#To save pickle file with full_output and pickle with full results
	pickle.dump(full_output, open(direc_name+'results_full_outputN'+str(N)+'.pkl', 'wb'))
	pickle.dump(FullResults, open(direc_name+'FullResultsN'+str(N)+'.pkl', 'wb'))

	#Find 1sig confidence interval and median
	median = np.median(all_depths, axis=0)
	lower = np.percentile(all_depths, 16, axis=0)
	upper = np.percentile(all_depths, 84, axis=0)

	# to add the offset to the right wavelength region of the retrieval models
	if InstNum:
		Cnt = 1
		for I in range(InstNum): 
			AddOff = "median[np.logical_and(wavLen > Params0["+str(InsCnt+Cnt+1)+"]," 
			AddOff += "wavLen < Params0["+str(InsCnt+Cnt+2)+"])] -= np.median(FullResults[:,"+str(InsCnt+Cnt)+"])" # - offset because in the fitting added the offset to the depth
			exec(AddOff) # for the best fit of the model
			AddOff_low = "lower[np.logical_and(wavLen > Params0["+str(InsCnt+Cnt+1)+"],"
			AddOff_low += "wavLen < Params0["+str(InsCnt+Cnt+2)+"])] -= np.median(FullResults[:,"+str(InsCnt+Cnt)+"])" 
			exec(AddOff_low) # for the 1sig lower error bar of the model
			AddOff_upper = "upper[np.logical_and(wavLen > Params0["+str(InsCnt+Cnt+1)+"],"
			AddOff_upper += "wavLen < Params0["+str(InsCnt+Cnt+2)+"])] -= np.median(FullResults[:,"+str(InsCnt+Cnt)+"])" 
			exec(AddOff_upper) # for the 1sig upper error bar of the
			Cnt +=3	

	#To convert to appropriate units (ppm Depth/RpRs)
	if PlotDepth: #*1e6 to convert to ppm
		lower = lower*1e6 
		median = median*1e6
		upper = upper*1e6
	else: # np.sqrt() to convert from depth to Rp/Rs
		lower = np.sqrt(lower) 
		median = np.sqrt(median)
		upper = np.sqrt(upper)
	
	# to convert from SI to microns
	wavLen=wavLen*1e6

	#to save best fit to .dat file
	Datafile = open(direc_name+'RetrievalFits'+SaveAs+'_N'+str(N)+'.dat', 'w')
	Datafile.write('#Wavelength		best fit	lower 1sig		upper 1sig	\n')
	for r in range(len(median)):
		Datafile.write(str(wavLen[r])+'	'+str(median[r])+'	'+str(lower[r])+'	'+str(upper[r])+'\n')
	Datafile.close()

whereWave = np.where(wavLen < MaxWav+.2)[0] # wavelength bounds for plotting


#######################To make and plot corner plots ########################
# To convert Rp to Rp/Rj, Mp/Mj, and Rs/Rsun from SI units
if 'Rp' in param_names:
	Rp_indx = param_names.index('Rp')
	results['samples'][:,Rp_indx] = results['samples'][:,Rp_indx]/R_jup

if 'Mp' in param_names:
	Mp_indx = param_names.index('Mp')
	results['samples'][:,Mp_indx] = results['samples'][:,Mp_indx]/M_jup

if 'Rs' in param_names:
	Rs_indx = param_names.index('Rs')
	results['samples'][:,Rs_indx] = results['samples'][:,Rs_indx]/R_sun

#To plot results in Corner plots
fig = corner.corner(results['samples'], weights=np.exp(results['logwt']), title_fmt=".5f",
                    range=[0.99] * results['samples'].shape[1], show_titles = True,
                    labels=param_names, quantiles = [0.15865,0.50, 0.84135])
fig.savefig(direc_name+"multinest_corner.png")
print ("saved corner plot in dir'"+direc_name+ "as 'multinest_corner.png'")
plt.close()

########################To plot best retrieval fit########################
#general settings
xlabel = 'Wavelength ($\mu$m)'
ylabel = 'Transit depth (ppm)'
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
rcParams['axes.linewidth'] = 1.2 
rcParams['xtick.direction'] = 'out'
rcParams['ytick.direction'] = 'out'
rcParams['lines.markeredgewidth'] = 1 
bin_kwargs = {"zorder":5}
model_kwargs = {"zorder":0}
error_kwargs = {"zorder":10}
fig = plt.figure(figsize=(19,10))#plt.figure(figsize=(25,15))
ax = fig.add_subplot(111)
facecols = ["white","#a31621","#4e8098","#faa916"]
edgecols = ["black","#a31621","#4e8098","#faa916"]

# Plot datasets:
for i in Insturments.keys():
	offset = 0
	instr_count = Insturments[i]['count']-1
	if InstNum:
		ParNum = InsCnt + ((instr_count*2)+1) # to get right offset based on the insturments order in .dat file, remember dics have no specific order 
		offset = -np.median(FullResults[:,ParNum])*1e6 # in depth. 1e6 to convert to ppm. Minus because retrieval added the offset when calculating depth. Thus, best depth is -offset
	w = ((Insturments[i]['wav up']+Insturments[i]['wav low'])/2.0)*1e-4 # convert from angst to microns
	w_errUP = (Insturments[i]['wav up']*1e-4)-w
	w_errLw = w-(Insturments[i]['wav low']*1e-4)
	if PlotDepth:
		plt.errorbar(w,Insturments[i]['Depth']+offset,xerr=[w_errLw,w_errUP],\
		    yerr=[Insturments[i]['Err dwn'],Insturments[i]['Err up']],\
		    fmt='o',markerfacecolor=facecols[instr_count],markeredgecolor=edgecols[instr_count],\
		    ecolor=edgecols[instr_count],elinewidth=1,markersize=7,label=i,**error_kwargs)
	else:
		RpRs = np.sqrt((Insturments[i]['Depth']+offset)*1e-6) #*1e-6 is to convert depth from ppm to just depth
		RpRs_errUp = np.sqrt((Insturments[i]['Depth']+offset+Insturments[i]['Err up'])*1e-6) - RpRs
		RpRs_errDwn = RpRs - np.sqrt((Insturments[i]['Depth']+offset-Insturments[i]['Err dwn'])*1e-6)
		plt.errorbar(w,RpRs,xerr=[w_errLw,w_errUP],yerr=[RpRs_errDwn,RpRs_errUp],fmt='o',\
		    markerfacecolor=facecols[instr_count],markeredgecolor=edgecols[instr_count],\
		    ecolor=edgecols[instr_count],elinewidth=1,markersize=7,label=i,**error_kwargs)
	instr_count += 1

#print Log Evidence in title
logZ = np.median(results['logz'])
logZerr = np.nanmedian(results['logzerr']) # the 1st value is nan for some reason???
plt.title(direc_name+'  (ln) Evidence: '+str(logZ)+'+/-'+str(logZerr))

#Plot 1sig confidence interval and best fit
if PlotDepth: #if you want to plot in terms of Depth (ppm) or RpRs
	yLabel = 'Transit depth (ppm)'
	yLim = [MinDepth-500, MaxDepth+500]
else:
	yLabel = 'Rp/Rs'
	yLim = [np.sqrt((MinDepth-500)*1e-6), np.sqrt((MaxDepth+500)*1e-6)]
plt.plot(wavLen[whereWave],median[whereWave], 'k', label = 'median retrieval model', linewidth=0.5)
plt.fill_between(wavLen[whereWave], lower[whereWave], upper[whereWave], color='cornflowerblue', alpha=0.2,edgecolor="none")
plt.ylim([yLim[0],yLim[1]])
plt.ylabel(yLabel)
plt.xlim([MinWav-.01, MaxWav+.1])
plt.xlabel('Wavelength ($\mu$m)')
plt.xscale('log')
plt.legend(ncol=3)

#To bin retrieval to resolution of data
binned_wavs = []
binned_depths = []
for (start, end) in bins:
	cond = np.where((wavLen >= start)&(wavLen < end))[0]
	binned_wavs.append(np.mean(wavLen[cond]))
	binned_depths.append(np.mean(median[cond]))
plt.plot(binned_wavs,binned_depths, 'ks')	

# Supress scientific notation:
from matplotlib.ticker import ScalarFormatter
for axis in [ax.xaxis]:
	axis.set_major_formatter(ScalarFormatter())
	axis.set_minor_formatter(ScalarFormatter())

plt.savefig(direc_name+'FinalTransSpect'+SaveAs+'_N'+str(N)+'.png')
print "Saved plot named 'FinalTransSpect"+SaveAs+"_N"+str(N)+"png' in '"+direc_name+"'\n"
ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
plt.tight_layout()
plt.show()
plt.close()

