from __future__ import print_function

import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate
import emcee
# import corner
import pickle
import retrieval_input as RI
import datetime
import time
from platon.fit_info import FitInfo
from platon.retriever import Retriever
from platon.constants import R_sun, R_jup, M_jup

if not os.path.exists(RI.Target): #To make the main directory of the retrieval results. Titled by the name of the target
    os.system('mkdir %s'%RI.Target)
direc_name = RI.Target+'/'

#Parameters form retrieval_input. take values of RI now because don't want to be reading in values from RI much later after starting script
SubFoldExtra = RI.SubFoldExtra
rp = RI.Rp
rpBnds = RI.RpBnds
mp = RI.Mp
mpErr = RI.MpErr
t_eq = RI.T_eq
t_eqBnds = RI.T_eqBnds
r_star = RI.R_star
r_starErr = RI.R_starErr
t_star = RI.T_star
TsErr = RI.TsErr
lg_scat_fact= RI.log_scatt_factor
lg_scat_factBnds = RI.log_scatt_factorBnds
lgZ = RI.logZ
lgZBnds = RI.logZBnds
lg_cloudtop_P = RI.log_cloudtop_P
lg_cloudtop_PBnds = RI.log_cloudtop_PBnds
CO_rat = RI.CO_ratio
CO_ratBnds = RI.CO_ratioBnds
scat_slope = RI.scatt_slope
scat_slopeBnds = RI.scatt_slopeBnds
err_multiple = RI.error_multiple
err_multipleBnds = RI.error_multipleBnds
t_spot = RI.T_spot
t_spotBnds = RI.T_spotBnds
spotCovFrac = RI.spot_cov_frac
spotCovFracBnds = RI.spot_cov_fracBnds
nlive = RI.nlive
NumInst = RI.NumInst
if NumInst:
    for o in range(NumInst): #To add the offset parameters -initial offset, offset uniform bounds, & wavelength range of offset insturment- to respective variables  
        exec("OffSet%s=%s" % (o+1,'RI.Offset'+str(o+1)))
        exec("OffSet%s_Rng=%s" % (o+1,'RI.Offset'+str(o+1)+'_Rng'))
        exec("OffSet%s_Wavs=%s" % (o+1,'RI.Offset'+str(o+1)+'_Wavs'))

# First, extract data:
data_fname = 'Dats/'+RI.dataset+'.dat'
print("Using: "+RI.dataset+'.dat')
data = np.genfromtxt(data_fname, dtype=None, encoding=None)

# Save it into our data dictionary:
data_dictionary = {}
instCnt = 1
for i in range(len(data)):
    wlow,wup,depth,errup,errlow,instrument,offset = data[i]
    if instrument not in data_dictionary.keys():
        data_dictionary[instrument] = {}
        data_dictionary[instrument]['count'] = instCnt # to keep track of the instrument order
        instCnt += 1
        if offset.lower() == 'yes' or offset.lower() == 'true':
            data_dictionary[instrument]['offset'] = True
        else:
            data_dictionary[instrument]['offset'] = False
        data_dictionary[instrument]['wlow'] = np.array([])
        data_dictionary[instrument]['wup'] = np.array([])
        data_dictionary[instrument]['depth'] = np.array([])
        data_dictionary[instrument]['depth_err'] = np.array([])
    data_dictionary[instrument]['wlow'] = np.append(data_dictionary[instrument]['wlow'],np.float(wlow))
    data_dictionary[instrument]['wup'] = np.append(data_dictionary[instrument]['wup'],np.float(wup))
    data_dictionary[instrument]['depth'] = np.append(data_dictionary[instrument]['depth'],np.float(depth))
    data_dictionary[instrument]['depth_err'] = np.append(data_dictionary[instrument]['depth_err'],(errup+errlow)/2.)
for InSt in data_dictionary.keys():
    data_dictionary[InSt]['MeanDepth'] = np.mean(data_dictionary[InSt]['depth'])

#To combine all the data together because platon doesn't care which insturment it came from
bins_up =  np.array([])
bins_dwn = np.array([])
depths = np.array([])
errors = np.array([])
Inst = list(data_dictionary.keys())
InstWavRng = []
for i in range(len(Inst)):
    bins_up = np.append(bins_up, data_dictionary[Inst[i]]['wup']*1e-10) # to convert from angstrom to meters
    bins_dwn = np.append(bins_dwn, data_dictionary[Inst[i]]['wlow']*1e-10) 
    depths = np.append(depths, data_dictionary[Inst[i]]['depth']*1e-6) # to convert from ppm to fraction
    errors = np.append(errors, data_dictionary[Inst[i]]['depth_err']*1e-6)
bins = np.zeros((len(bins_up),2))
bins[:,0] = bins_dwn
bins[:,1] = bins_up

#create a Retriever object
retriever = Retriever()

#create a FitInfo object and set best guess parameters
FitInfoStr = "Rs=r_star*R_sun, Mp=mp*M_jup, Rp=rp*R_jup, T=t_eq,logZ=lgZ, CO_ratio=CO_rat, T_star=t_star,"
FitInfoStr += "spot_cov_frac=spotCovFrac, T_spot=t_spot,log_cloudtop_P =lg_cloudtop_P, InstNum = NumInst," 
if NumInst:
    for n in range(NumInst): # to add the N number of offsets
        FitInfoStr += 'Offset'+str(n+1)+'=OffSet'+str(n+1)+','+'Offset'+str(n+1)+'_Wavs=OffSet'+str(n+1)+'_Wavs,' 
exec("fit_info = retriever.get_default_fit_info(%s)" % (FitInfoStr))

#Will create a .txt file with the details of the fitting
txtInfo = "This specific retrieval model was done by:\n"

#Add fitting parameters - this specifies which parameters you want to fit
direc_name += "FIT_"
if r_starErr:
    fit_info.add_gaussian_fit_param('Rs', r_starErr*R_sun)
    direc_name += "Rs"
    txtInfo += "*Fitting radius of the star (Rs) with gaussian priors containing a standard deviation of "+str(r_starErr)+"R_sun.\n"
if mpErr:
    fit_info.add_gaussian_fit_param('Mp', mpErr*M_jup)
    direc_name += "Mp"
    txtInfo += "*Fitting mass of the planet (Mp) with gaussian priors containing a standard deviation of "+str(mpErr)+"M_jup.\n"
if rpBnds:
    fit_info.add_uniform_fit_param('Rp', rpBnds[0]*R_jup, rpBnds[1]*R_jup)
    direc_name += "Rp"
    txtInfo += "*Fitting radius of the planet (Rp) with uniform priors containing bounds of "+str(rpBnds[0])+"R_jup to "+str(rpBnds[1])+"R_jup.\n"
if t_eqBnds:
    fit_info.add_uniform_fit_param('T', t_eqBnds[0], t_eqBnds[1])
    direc_name += "Tp"
    txtInfo += "*Fitting equilibrium Temp of the star (Tp) with uniform priors containing bounds of "+str(t_eqBnds[0])+"K to "+str(t_eqBnds[1])+"K.\n"
if lgZBnds:
    fit_info.add_uniform_fit_param("logZ", lgZBnds[0], lgZBnds[1])
    direc_name += "Lz"
    txtInfo += "*Fitting log metallicity of the star (Lz) with uniform priors containing bounds of "+str(lgZBnds[0])+"dex to "+str(lgZBnds[1])+"dex.\n"
if CO_ratBnds:
    fit_info.add_uniform_fit_param("CO_ratio",CO_ratBnds[0],CO_ratBnds[1])
    direc_name += "Co"
    txtInfo += "*Fitting CO ratio of the planet (Co) with uniform priors containing bounds of "+str(CO_ratBnds[0])+" to "+str(CO_ratBnds[1])+".\n"
if lg_cloudtop_PBnds:
    fit_info.add_uniform_fit_param("log_cloudtop_P", lg_cloudtop_PBnds[0], lg_cloudtop_PBnds[1])
    direc_name += "Lp"
    txtInfo += "*Fitting log cloud top Pressure (Lp) with uniform priors containing bounds of "+str(lg_cloudtop_PBnds[0])+"dex to "+str(lg_cloudtop_PBnds[1])+"dex.\n"
if lg_scat_factBnds:
    fit_info.add_uniform_fit_param("log_scatt_factor", lg_scat_factBnds[0], lg_scat_factBnds[1])
    direc_name += "Ls"
    txtInfo += "*Fitting log scattering factor (Ls) with uniform priors containing bounds of "+str(lg_scat_factBnds[0])+"dex to "+str(lg_scat_factBnds[1])+"dex.\n"
if scat_slopeBnds:
    fit_info.add_uniform_fit_param("scatt_slope",scat_slopeBnds[0],scat_slopeBnds[1])
    direc_name += "Ss"
    txtInfo += "*Fitting scattering slope (Ss) with uniform priors containing bounds of "+str(scat_slopeBnds[0])+" to "+str(scat_slopeBnds[1])+".\n"
if err_multipleBnds:
    fit_info.add_uniform_fit_param("error_multiple", err_multipleBnds[0], err_multipleBnds[1])
    direc_name += "Em"
    txtInfo += "*Fitting error multiple (Em) with uniform priors containing bounds of "+str(err_multipleBnds[0])+" to "+str(err_multipleBnds[1])+".\n"   
if TsErr:
    fit_info.add_gaussian_fit_param('T_star', TsErr)
    direc_name += "Ts"
    txtInfo += "*Fitting Temp of the star (Ts) with gaussian priors containing a standard deviation of "+str(TsErr)+"K.\n"
if t_spotBnds:
    fit_info.add_uniform_fit_param("T_spot", t_spotBnds[0],t_spotBnds[1])
    direc_name += "Th"
    txtInfo += "*Fitting Temp of stellar spots/heterogeneities (Th) with uniform priors containing bounds of "+str(t_spotBnds[0])+"K to "+str(t_spotBnds[1])+"K.\n"   
if spotCovFracBnds:
    fit_info.add_uniform_fit_param("spot_cov_frac", spotCovFracBnds[0], spotCovFracBnds[1])
    direc_name += "Cf"
    txtInfo += "*Fitting stellar spots covering fraction (Cf) with uniform priors containing bounds of "+str(spotCovFracBnds[0])+" to "+str(spotCovFracBnds[1])+".\n"   
if NumInst:
    direc_name += 'Os'   
    for I in data_dictionary.keys():
        if data_dictionary[I]['offset']:
            exec("BndOff = data_dictionary[I]['MeanDepth']*OffSet"+str(data_dictionary[I]['count'])+"_Rng*1e-6")
            AddUnifStr = "'"+'Offset'+str(data_dictionary[I]['count'])+"',-BndOff, BndOff"
            exec("fit_info.add_uniform_fit_param(%s)"%(AddUnifStr))
            txtInfo += "*Fitting depth offset for "+I+" data with uniform priors containing bounds of "+str(-BndOff)+" to "+str(BndOff)+".\n"
if nlive != 1000: # the default nlives is 1000
    direc_name += 'N'+str(nlive)
txtInfo += "The number of live points (number of samples drawn per layer) was " +str(nlive)+".\n"  #want to always print number of live points used
if SubFoldExtra: #To add extra string values to general subfolder naming scheme
    direc_name += SubFoldExtra+'/'
else:
    direc_name += "/"

# To make results output directories
if not os.path.exists(direc_name): 
    os.system('mkdir %s'%direc_name)
os.system("cp %s %s"%(data_fname, direc_name))

#To print summary .txt file
now = datetime.datetime.now()
StartT = time.time()
TimeStamp = str(now.year)+'-'+str(now.month)+'-'+str(now.day)
SummaryF = open(direc_name+'INFO'+TimeStamp+'.txt','w+')
SummaryF.write("Run started on "+TimeStamp+' '+str(now.hour)+'.'+str(now.minute)+'.'+str(now.second)+'\n\n\n')
SummaryF.write(txtInfo)
SummaryF.close()

#Use Nested Sampling to do the fitting
result = retriever.run_multinest(bins, depths, errors, fit_info, plot_best=False,npoints=nlive) 
pickle.dump(result,open(direc_name+'multinest_result.pkl','wb'))
pickle.dump(fit_info,open(direc_name+'multinest_FitInfo.pkl','wb'))

#To output final timestamps in summary .txt file
End = datetime.datetime.now()
TotalT = (time.time()-StartT)/(3600.0) #Total run time in hours
SummaryF = open(direc_name+'INFO'+TimeStamp+'.txt','a+')
SummaryF.write('\n\n\n')
SummaryF.write("Run ended on "+str(End.year)+'-'+str(End.month)+'-'+str(End.day)+' '+str(End.hour)+'.'+str(End.minute)+'.'+str(End.second)+'\n')
SummaryF.write("Total run time: "+str(TotalT)+"hr \n")
SummaryF.close()

