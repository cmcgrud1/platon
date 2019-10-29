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

def write_param_estimates_file(samples, best_params, best_lnprob, fit_labels,filename="BestFit.txt"):
	#To save best fits of retreival. Same function in Mike Zhang's _output_writer.py script
	output = "#Parameter Lower_error Median Upper_error Best_fit\n"
	output += "Max_lnprob {}\n".format(best_lnprob)
	BESTFIT = {} #To bestfit parameters
	for i, name in enumerate(fit_labels):
		lower = np.percentile(samples[:, i], 16)
		median = np.median(samples[:, i])
		upper = np.percentile(samples[:, i], 84)
		best = best_params[i]
		output += "{} {} {} {} {}\n".format(name, \
			median - lower, median, upper - median, best)
		BESTFIT[name] = np.array([median-lower, median, upper-median, best])
	print(output)
	with open(filename, "w") as f:
		f.write(output)
	return BESTFIT

#To save plot settings
N = PRI.N # number of elements for each chain to resample, advised 1000
if N < 10:
	N = 10 # can't have less than 10 samples
PlotDepth = PRI.PlotDepth #if you want to plot in terms of Depth (ppm) or RpRs
X_bounds = PRI.X_bounds
Y_bounds = PRI.Y_bounds

#To load info that the retriever code used
direc_name = PRI.Target+'/'+PRI.SubFolder+'/'

#List of all compute_depths input parameters, in the right order
ParamNames = ['Rs', 'Mp', 'Rp', 'T', 'logZ', 'CO_ratio', 'log_cloudtop_P', 'log_scatt_factor', 'scatt_slope', 'T_star', 'T_spot', 'spot_cov_frac'] # error_multiple can be fit for in ln_like, but that's after computing depth 
InstNum = PRI.InstNum
Params0 =[PRI.R_star*R_sun, PRI.Mp*M_jup, PRI.Rp*R_jup, PRI.T_eq, PRI.logZ, PRI.CO_ratio, PRI.log_cloudtop_P, PRI.log_scatt_factor, PRI.scatt_slope, PRI.T_star, PRI.T_spot, PRI.spot_cov_frac]

#To add offset terms if they exist
if InstNum:
	ParamNames.append('InstNum')
	Params0.append(InstNum)
	for i in range(InstNum):
		ParamNames.append('Offset'+str(i+1))
		exec("Params0.append(PRI.Offset"+str(i+1)+")")
		ParamNames.append('Offset'+str(i+1)+'Dwn_Wavs') # Need to seperate upper and lower wavelengths of insturment because need to easily put values in array
		exec("Params0.append(PRI.Offset"+str(i+1)+"_Wavs[0])")
		ParamNames.append('Offset'+str(i+1)+'Up_Wavs') #Not same name as in compute_depths, but doesn't matter because will always be fixed anyways
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
Results1Exist = os.path.isfile(direc_name+'bestfit_full_output.pkl') #This is the number of samples for resampling DIFFERENT FROM NLIVE!
Results2Exist = os.path.isfile(direc_name+'HigResRetrievTransSpec'+SaveAs+'_N'+str(N)+'.dat')
Results3Exist = os.path.isfile(direc_name+'DataResRetrievTransSpec'+SaveAs+'_N'+str(N)+'.dat')
Results4Exist = os.path.isfile(direc_name+'BestFit.txt')
if Results1Exist and Results2Exist and Results3Exist and Results4Exist:
	print ("Found 'bestfit_full_output.pkl', '<Hig/Data>ResRetrievTransSpec"+SaveAs+"_N"+str(N)+".dat', and BestFit.txt file.")
	print ("Using '"+direc_name+"<Hig/Data>ResRetrievTransSpecs"+SaveAs+"_N"+str(N)+".dat' files to make plots... \n")
	wavLen,best, lower, upper = np.loadtxt(direc_name+'HigResRetrievTransSpec'+SaveAs+'_N'+str(N)+'.dat', unpack = True)
	binned_wavs, binned_depths = np.loadtxt(direc_name+'DataResRetrievTransSpec'+SaveAs+'_N'+str(N)+'.dat', unpack = True)
	BESTFIT = {}
	bestFits = open(direc_name+'BestFit.txt', 'r')
	for line in bestFits:
		if line.startswith('#'):
			pass
		else:
			split = line.split() 
			if len(split) > 2: # don't want max lnProb, which only has 2 elements in that line
				BESTFIT[split[0]] = np.array([np.float(split[1]), np.float(split[2]), np.float(split[3]), np.float(split[4])])

#Only do resampling if hasn't been done for this set up
else:	
	print ("Didn't find 'bestfit_full_outputN"+str(N)+".pkl', '<Hig/Data>ResRetrievTransSpec"+SaveAs+"_N"+str(N)+".dat', and BestFit.txt file.")
	print ("Resampling transit depth parameters based on '"+direc_name+"multinest_result.pkl' file... \n")
	#To pull out the sampled chains for each fitted for parameter
	chain = dyutil.resample_equal(results['samples'], results['weights'])

	#To get parameters in right arrangement needed for compute_depths(). i.e put all the fixed params in
	#For the values that were not fitted for set equal to the initial values (Params0)
	FullResults = np.ones((results['samples'].shape[0],len(ParamNames)))
	print ("results['samples'].shape", results['samples'].shape)
	print (' param_names:',  param_names)
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
	Rand_N = list(np.random.choice(FullResults.shape[0],N, replace = False))
	for params in FullResults[Rand_N][:]: #To select N random chains
		if (step%prntStp) == 0: # to print update
			print ('step:'+str(step)+'/'+str(N-1))
		cld_tp_P = 10**params[6] #covert log_cloudtop_press to cloudtop_press
		sct_fct = 10**params[7] #covert log_scatt_fract to scatt_fract		
		CmpDept_kwargs = {"logZ":params[4], "CO_ratio":params[5], "cloudtop_pressure":cld_tp_P, "scattering_factor":sct_fct,\
		"scattering_slope":params[8], "T_star":params[9], "T_spot":params[10], "spot_cov_frac":params[11]}
		if InstNum: #to add offsets if need be
			CmpDept_kwargs["InstNum"] = InstNum
			cnt =1
			for I_N in range(InstNum):
				CmpDept_kwargs['Offset'+str(I_N+1)] = params[InsCnt+cnt]
				CmpDept_kwargs['Offset'+str(I_N+1)+'_Wavs']=[params[(InsCnt+cnt+1)],params[(InsCnt+cnt+2)]] 
				cnt += 3
		# To compute the transit spectrum and full_output parameters for each chain
		wavLen, depths = calculator.compute_depths(params[0],params[1],params[2],params[3], **CmpDept_kwargs) 
		all_depths[step] = depths #To save all the depth values of the resampling
		step += 1
	all_depths = np.array(all_depths)

	#To create 'Best_fit.txt' file. It's filled with the Max lnProb and median, 1sig confidence interval, and best fit of all fitted parameters
	bst_parm_ar = results['samples'][np.argmax(results['logp'])]
	print ('\n')
	BESTFIT = write_param_estimates_file(chain, bst_parm_ar, np.max(results['logp']), param_names, filename=direc_name+"BestFit.txt")
	
	#To calculate the full_output for the best fit parameters. Had to use the best fit parameters to calculate the full_output 
	#because the shape of the full_output arrays will be different based on what cloud top pressure is used, which changes.
	#Thus, can't take median or confidence bars of that
	BstParms = Params0
	for bf in BESTFIT.keys():
		if bf == 'error_multiple': #still don't know how to properly include error_multiple term
			pass
		else:
			BstParms[ParamNames.index(bf)] = BESTFIT[bf][3] #To get the best fit parameters
	cld_tp_P, sct_fct= 10**BstParms[6],10**BstParms[7]  #covert from log
	CmpDept_kwargs = {"logZ":BstParms[4], "CO_ratio":BstParms[5], "cloudtop_pressure":cld_tp_P, "scattering_factor":sct_fct,\
	"scattering_slope":BstParms[8], "T_star":BstParms[9], "T_spot":BstParms[10], "spot_cov_frac":BstParms[11],"full_output":True}
	if InstNum: #to add offsets if need be
		CmpDept_kwargs["InstNum"] = InstNum
		cnt =1
		for I_N in range(InstNum):
			CmpDept_kwargs['Offset'+str(I_N+1)] = BstParms[InsCnt+cnt]
			CmpDept_kwargs['Offset'+str(I_N+1)+'_Wavs']=[BstParms[(InsCnt+cnt+1)],BstParms[(InsCnt+cnt+2)]] 
			cnt += 3
	wavLen, best, full_output = calculator.compute_depths(BstParms[0],BstParms[1],BstParms[2],BstParms[3], **CmpDept_kwargs) 

	#To save pickle file with full_output and pickle with full results
	print ('dumping full_output...')
	pickle.dump(full_output, open(direc_name+'bestfit_full_output.pkl', 'wb'))

	#Find 1sig confidence intervals of transit depths. Technically this only works if the distribution is gaussian. However, we'll aprox. ;)
	# median = np.median(all_depths, axis=0) #using best fit instead of median. Because median is only best fit when distribution is gaussian. That's NEVER true
	lower = np.percentile(all_depths, 16, axis=0)
	upper = np.percentile(all_depths, 84, axis=0)

	# to add the offset to the right wavelength region of the retrieval models
	if InstNum:
		Cnt = 1
		for I in range(InstNum):
			if 'Offset'+str(I+1) in BESTFIT:
				# for the best fit of the model
				cond1, cond2 = wavLen > Params0[InsCnt+Cnt+1], wavLen < Params0[InsCnt+Cnt+2]
				best[np.logical_and(cond1,cond2)] -= BESTFIT['Offset'+str(I+1)][3]# - offset because in the fitting added the offset to the depth				
				lower[np.logical_and(cond1,cond2)] -= BESTFIT['Offset'+str(I+1)][3]# for the 1sig lower error bar of the model
				upper[np.logical_and(cond1,cond2)] -= BESTFIT['Offset'+str(I+1)][3] # for the 1sig upper error bar of the
			Cnt +=3	

	#To convert to appropriate units (ppm Depth or RpRs)
	if PlotDepth: #*1e6 to convert to ppm
		lower = lower*1e6 
		best = best*1e6
		upper = upper*1e6
	else: # np.sqrt() to convert from depth to Rp/Rs
		lower = np.sqrt(lower) 
		best = np.sqrt(best)
		upper = np.sqrt(upper)
	
	# to convert from SI to microns
	wavLen=wavLen*1e6

	#To bin retrieval to resolution of data
	binned_wavs = []
	binned_depths = []
	for (start, end) in bins:
		cond = np.where((wavLen >= start)&(wavLen < end))[0]
		binned_wavs.append(np.mean(wavLen[cond]))
		binned_depths.append(np.mean(best[cond]))

	#to save best fit and 1sig confidence bars at high resolution to .dat file
	print ('Saved retrieval best fit transmission spectrum .dat files...')
	Datafile = open(direc_name+'HigResRetrievTransSpec'+SaveAs+'_N'+str(N)+'.dat', 'w')
	Datafile.write('#Wavelength		best fit	lower 1sig		upper 1sig	\n')
	for r in range(len(best)):
		Datafile.write(str(wavLen[r])+'	'+str(best[r])+'	'+str(lower[r])+'	'+str(upper[r])+'\n')
	Datafile.close()

	#To save best fit at the data's resolution to .dat file
	Datafile = open(direc_name+'DataResRetrievTransSpec'+SaveAs+'_N'+str(N)+'.dat', 'w')
	Datafile.write('#Wavelength		best fit \n')
	for b in range(len(binned_depths)):
		Datafile.write(str(binned_wavs[b])+'	'+str(binned_depths[b])+'\n')
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
#TODO: have corner plot a vertical line that represents the best fit (different from the mean), with a different color from the confidence intervals
fig = corner.corner(results['samples'], weights=np.exp(results['logwt']), title_fmt=".5f",
                    range=[0.99] * results['samples'].shape[1], show_titles = True,
                    labels=param_names, quantiles = [0.15865,0.50, 0.84135])
fig.savefig(direc_name+"multinest_corner.png")
print ("Saved corner plot in dir '"+direc_name+ "as 'multinest_corner.png'")
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
for i in list(Insturments.keys()):
	offset = 0
	instr_count = Insturments[i]['count']-1
	if 'Offset'+str(instr_count+1) in BESTFIT:
		offset = -BESTFIT['Offset'+str(instr_count+1)][3]*1e6 # in depth. 1e6 to convert to ppm. Minus because retrieval added the offset when calculating depth. Thus, best depth is -offset
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
logZ = results['logz'][np.argmax(results['logp'])] #NOT SURE IF THIS IS RIGHT WAY TO GET LOGZ!!!
logZerr = results['logzerr'][np.argmax(results['logp'])] #NOT SURE IF THIS IS RIGHT WAY TO GET LOGZ err!!!
plt.title(direc_name+'  (ln) Evidence: '+str(logZ)+'+/-'+str(logZerr))

#Plot 1sig confidence interval and best fit
if PlotDepth: #if you want to plot in terms of Depth (ppm) or RpRs
	yLabel = 'Transit depth (ppm)'
	yLim = [MinDepth-500, MaxDepth+500]
else:
	yLabel = 'Rp/Rs'
	yLim = [np.sqrt((MinDepth-500)*1e-6), np.sqrt((MaxDepth+500)*1e-6)]
xLim = [MinWav-.01, MaxWav+.1]
if Y_bounds: #if you specified the plotting bounds you want
	yLim = Y_bounds
if X_bounds:
	xLim = X_bounds
plt.plot(binned_wavs,binned_depths, 'ks') #plot best fit at data resolution
plt.plot(wavLen[whereWave],best[whereWave], 'k', label = 'best retrieval model', linewidth=0.5)
plt.fill_between(wavLen[whereWave], lower[whereWave], upper[whereWave], color='cornflowerblue', alpha=0.2,edgecolor="none")
plt.ylim([yLim[0],yLim[1]])
plt.ylabel(yLabel)
plt.xlim([xLim[0],xLim[1]])
plt.xlabel('Wavelength ($\mu$m)')
plt.xscale('log')
plt.legend(ncol=3)

# Supress scientific notation:
from matplotlib.ticker import ScalarFormatter
for axis in [ax.xaxis]:
	axis.set_major_formatter(ScalarFormatter())
	axis.set_minor_formatter(ScalarFormatter())

plt.savefig(direc_name+'FinalTransSpect'+SaveAs+'_N'+str(N)+'.png')
print ("Saved plot named 'FinalTransSpect"+SaveAs+"_N"+str(N)+"png' in '"+direc_name+"'\n")
ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
plt.tight_layout()
plt.show()
plt.close()

