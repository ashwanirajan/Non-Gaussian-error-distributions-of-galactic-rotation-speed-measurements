import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.special import binom
from scipy import stats
import math
from decimal import *

#############importing data
filename = "galactic_mass_modelling.csv"
D_0 = np.loadtxt(filename, delimiter = ",", usecols= [2])
error_0 = np.loadtxt(filename, delimiter= ",", usecols= [3])
R_0 = np.loadtxt(filename, delimiter= ",", usecols= [4])
D_i = D_0*8.3/R_0
#D_i[22] = 212
#D_i[38] = 209
#D_i[76] = 214
error = error_0*8.3/R_0
#error[22] = 16
#error[38] = 2
#error[76] = 8


##############end importing


##############Defining Sigma

sigma_i = np.array(np.sqrt(error*error))

#print(sigma_i)
weights_list = []

for j in range(0,14):
		weights_list.append(1/(sigma_i[j])**2)
weights = np.array(weights_list)
#print(weights)

##############Sigma defined

##############weighted mean
D_wm = sum(D_i*weights)/sum(weights)

print("weighted mean" + ' = ' + str(D_wm))
sigma_wm = 1/np.sqrt(sum(weights))
print("standard deviation of errors" + ' = ' + str(sigma_wm))
##############

##############weighted median (average of the two middle terms in the sorted array)

ordered_D_i = np.sort(D_i)
ordered_sigma_i = np.sort(sigma_i)
#print(ordered_sigma_i)
#print(ordered_D_i)
D_med = 0.5*(ordered_D_i[6] + ordered_D_i[7]) 
print("median" + "=" + str(D_med))

##############finding error in median

for k in range(0, 14):
			sum=0
			for j in range(0, k):
				sum += binom(14, j)
				if sum >=(2**13*0.317) :
					break
#print(j)
sigma_med = 0.5*(ordered_D_i[15-j] - ordered_D_i[j-1])
print("median error " + "= " + str(sigma_med))

################


#################Defining N_sigma_i array

N_Sigma_i_wm = (D_i - D_wm)/ np.sqrt((sigma_i)**2 + (sigma_wm)**2) #N_sigma_i array
N_Sigma_i_med = (D_i - D_med)/ np.sqrt((sigma_i)**2 + (sigma_med)**2)




D_value_gaussian_wm, p_value_gaussian_wm = stats.kstest(N_Sigma_i_wm,'norm')
print(stats.kstest(N_Sigma_i_wm,'norm'))

print("for Gaussian, p value(wm) is " + str(p_value_gaussian_wm))


D_value_gaussian_med, p_value_gaussian_med = stats.kstest(N_Sigma_i_med,'norm')
print(stats.kstest(N_Sigma_i_med,'norm'))

print("for Gaussian, p value(med) is " + str(p_value_gaussian_med))