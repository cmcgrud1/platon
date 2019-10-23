import numpy as np
from astropy import constants

bf = np.genfromtxt('BestFit.txt',skip_header=2,usecols=[1,2,3,4]) # values
bf_k = np.genfromtxt('BestFit.txt',skip_header=2,usecols=0,dtype=str) # keys

for k,v in zip(bf_k,bf):

    v_low = v[0]
    v_med = v[1]
    v_up = v[2]
    best_fitting_value = v[3]

    if "T" in k: # print Teq as integer
        print("%s = $ %d ^{+%d} _{-%d} $ ; best fitting value = %d"%(np.round(k),np.round(v_med),np.round(v_up),np.round(v_low),np.round(best_fitting_value)))

    elif "Rp" in k: # print Rp in R_jup to 2 decimal places
        print("%s = $ %.2f ^{+%.2f} _{-%.2f} $ ; best fitting value = %.2f"%(np.round(k,2),np.round(v_med/constants.R_jup.value,2),\
        np.round(v_up/constants.R_jup.value,2),np.round(v_low/constants.R_jup.value,2),np.round(best_fitting_value/constants.R_jup.value,2)))

    else: # print all other variables to 2 decimal places
        print("%s = $ %.2f ^{+%.2f} _{-%.2f} $ ; best fitting value = %.2f"%(np.round(k,2),np.round(v_med,2),np.round(v_up,2),np.round(v_low,2),np.round(best_fitting_value,2)))

max_lnprob = np.genfromtxt('BestFit.txt',skip_header=1,skip_footer=len(bf))[1].astype(float)

print("-------")
print("Max lnprob = %.2f"%(np.round(max_lnprob,2)))
