import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.special import binom
from scipy import stats
import math
from decimal import *

#############importing data
filename = "rotation1.csv"
D_0 = np.loadtxt(filename, delimiter = ",", usecols= [2])
error_0 = np.loadtxt(filename, delimiter= ",", usecols= [3])
R_0 = np.loadtxt(filename, delimiter= ",", usecols= [4])
D_i = D_0*8.3/R_0
D_i[22] = 212
D_i[38] = 209
D_i[76] = 214
error = error_0*8.3/R_0
error[22] = 16
error[38] = 2
error[76] = 8


#print(D_0)
#print(D_i)
#print(error)
#print(R_0)

##############end importing


##############Defining Sigma

sigma_i = np.array(np.sqrt(error*error))

#print(sigma_i)
weights_list = []

for j in range(0,137):
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
D_med = ordered_D_i[68] 
print("median" + "=" + str(D_med))

##############finding error in median

for k in range(0, 137):
			sum=0
			for j in range(0, k):
				sum += binom(137, j)
				if sum >=((2**136)*0.317) :
					break
print(j)
sigma_med = 0.5*(ordered_D_i[138-j] - ordered_D_i[j-1])
print("median error " + "= " + str(sigma_med))

################

################ These 4 lines below choose the central estimate for analysis, for making mean as the central estimate activate the first two lines below.
#D_ce = D_wm
#sigma_ce = sigma_wm
#D_ce = D_med
#sigma_ce = sigma_med

################goodness of fit

#X = sum((D_i - D_ce)*(D_i - D_ce)*weights)/ 134-1
#print("goodness of fit = " + str(X))


################No. of standard deviations

#N_sigma = abs(np.sqrt(X)-1)*np.sqrt(2*(134-1))
#print("N_sigma = " + str(N_sigma))
#####################

#################Defining N_sigma_i array

N_Sigma_i_wm = (D_i - D_wm)/ np.sqrt((sigma_i)**2 + (sigma_wm)**2) #N_sigma_i array
N_Sigma_i_med = (D_i - D_med)/ np.sqrt((sigma_i)**2 + (sigma_med)**2)
#print(N_Sigma_i_med)
#print(np.amin(N_Sigma_i))
#print(np.amax(N_Sigma_i))

################

################Histogram 1
bins_wm=np.arange(np.amin(N_Sigma_i_wm) - 0.5,np.amax(N_Sigma_i_wm) + 0.5 , 0.5)
bins2_wm=np.arange(np.amin(N_Sigma_i_wm) - 0.5,np.amax(N_Sigma_i_wm) + 0.5 , 2)
#print(bins)
plt.hist(N_Sigma_i_wm,bins=bins_wm,alpha=0.4,histtype='stepfilled', normed = True, edgecolor = 'black', linewidth=0.8)
plt.xlabel("$N_σ$")
plt.ylabel("probability density")

##############fitting a gaussian to the above histogram
parameters = norm.fit(N_Sigma_i_wm)
pdf_x_wm = np.linspace(np.amin(N_Sigma_i_wm) - 0.5,np.amax(N_Sigma_i_wm) + 0.5 ,500)
fitted_pdf_wm = norm.pdf(pdf_x_wm,loc = parameters[0],scale = parameters[1])

plt.plot(pdf_x_wm,fitted_pdf_wm,"black", linestyle="dashed", linewidth=1.5)
plt.legend()

plt.show()

###############histogram for median
bins_med=np.arange(np.amin(N_Sigma_i_med) - 0.5,np.amax(N_Sigma_i_med) + 0.5 , 0.5)
bins2_med=np.arange(np.amin(N_Sigma_i_med) - 0.5,np.amax(N_Sigma_i_med) + 0.5 , 2)
#print(bins)
binned_array, b ,c = plt.hist(N_Sigma_i_med,bins=bins_med,alpha=0.4,histtype='stepfilled', normed = True, edgecolor = 'black', linewidth=0.8)
plt.xlabel("$N_σ$")
plt.ylabel("probability density")

###################




parameters = norm.fit(N_Sigma_i_med)
pdf_x_med = np.linspace(np.amin(N_Sigma_i_med) - 0.5,np.amax(N_Sigma_i_med) + 0.5 ,500)
fitted_pdf_med = norm.pdf(pdf_x_med,loc = parameters[0],scale = parameters[1])

plt.plot(pdf_x_med,fitted_pdf_med,"black", linestyle="dashed", linewidth=1.5)
plt.legend()
plt.show()
###############

ordered_N_sigma_wm = np.sort(N_Sigma_i_wm)
ordered_N_sigma_med = np.sort(N_Sigma_i_med)
print(ordered_N_sigma_wm)
print(ordered_N_sigma_med)




plt.show()







###############Symmetric absolute value histogram 

abs_N_sigma_wm = np.abs(N_Sigma_i_wm)
N_sigma_symm_wm = []


b =  len(abs_N_sigma_wm)

###############this is for creating the N_sigma dataset for symmetric histogram
for u in range(0, b)	:
		if abs_N_sigma_wm[u] == -0. :
			N_sigma_symm_wm.append(-0.00001)
		else	:
			N_sigma_symm_wm.append(-1*abs_N_sigma_wm[u])

for v in range(0, b)	:
		N_sigma_symm_wm.append(abs_N_sigma_wm[v])

N_sigma_symm_arr_wm= np.array(N_sigma_symm_wm)
#print(N_sigma_symm_arr)

plt.hist(N_sigma_symm_arr_wm,bins=np.arange(-10-0.5,10+0.5,0.5),alpha=0.4,histtype='stepfilled', normed = True, weights=0.5*np.ones(len(N_sigma_symm_arr_wm)), edgecolor = 'black', linewidth=0.8)
plt.xlabel("$|N_σ|$")
plt.ylabel("probability density")


##############fitting a gaussian to the symmetric absolute histogram
parameters1 = norm.fit(N_sigma_symm_arr_wm)
pdf_x1_wm = np.linspace(-10 - 0.5,10 + 0.5 ,500)
fitted_pdf1_wm = norm.pdf(pdf_x1_wm,loc = parameters1[0],scale = parameters1[1])
plt.plot(pdf_x1_wm,fitted_pdf1_wm,"black",linestyle="dashed", linewidth=1.5)
plt.legend()


plt.show()
###############



###############Symmetric absolute value histogram for median
abs_N_sigma_med = np.abs(N_Sigma_i_med)

N_sigma_symm_med = []


b =  len(abs_N_sigma_med)

###############this is for creating the N_sigma dataset for symmetric histogram
for u in range(0, b)	:
		if abs_N_sigma_med[u] == -0. :
			N_sigma_symm_med.append(-0.00001)
		else	:
			N_sigma_symm_med.append(-1*abs_N_sigma_med[u])

for v in range(0, b)	:
		N_sigma_symm_med.append(abs_N_sigma_med[v])

N_sigma_symm_arr_med= np.array(N_sigma_symm_med)
#print(N_sigma_symm_arr)

plt.hist(N_sigma_symm_arr_med,bins=np.arange(-10-0.5,10+0.5,0.5),alpha=0.4,histtype='stepfilled', normed = True, weights=0.5*np.ones(len(N_sigma_symm_arr_med)), edgecolor = 'black', linewidth=0.8)
plt.xlabel("$|N_σ|$")
plt.ylabel("probability density")


##############fitting a gaussian to the symmetric absolute histogram
parameters1 = norm.fit(N_sigma_symm_arr_med)
pdf_x1_med = np.linspace(-10 - 0.5,10 + 0.5 ,500)
fitted_pdf1_med = norm.pdf(pdf_x1_med,loc = parameters1[0],scale = parameters1[1])
plt.plot(pdf_x1_med,fitted_pdf1_med,"black",linestyle="dashed", linewidth=1.5)
plt.legend()
###############







plt.show()









#####################for figure 2(comparison b/w mean median and expected gaussian)

N_Sigma_i1 = (D_i - D_wm)/ np.sqrt((sigma_i)**2 + (sigma_wm)**2)
plt.hist(np.abs(N_Sigma_i1),bins=np.arange(0, 10, 0.1),alpha=1,histtype='step', normed = True, edgecolor = 'blue', linestyle="-.", linewidth=1.7, label= "weighted mean")

N_Sigma_i2 = (D_i - D_med)/ np.sqrt((sigma_i)**2 + (sigma_med)**2)
plt.hist(np.abs(N_Sigma_i2),bins=np.arange(0, 10, 0.1),alpha=1,histtype='step', normed = True, edgecolor = 'red', linestyle="--", linewidth=1.7, label = "median")
plt.hist(np.random.randn(10000000),bins=np.arange(0, 10, 0.1),alpha=1,histtype='step', normed = True, edgecolor = 'black', linestyle="solid", linewidth=1.2, label = "expected value")
plt.xlabel("$|N_σ|$", fontsize = 15)
plt.ylabel("No. of measurements", fontsize = 16)
plt.legend(fontsize = 14, labelspacing = 0.9)
plt.show()
###############







plt.show()
 
######ks test input binned case 1. mean
binned_array_wm, b ,c = plt.hist(N_Sigma_i_wm,bins=bins_wm,alpha=0.4,histtype='stepfilled', normed = False, edgecolor = 'black', linewidth=0.8)
#print(binned_array_wm)
bin_no = len(binned_array_wm)
print(bin_no)
print(bins_wm)
central_value_wm = []

for s in range(1, bin_no+1) :
	central_value_wm.append(0.5*(bins_wm[s-1] + bins_wm[s]))
central_value_arr_wm = np.array(central_value_wm)
#print(central_value_arr_wm)

bin_input_wm = []
for h in range(0, 43) :
	for i in range(0, int(binned_array_wm[h])) :
		bin_input_wm.append(central_value_arr_wm[h])

bin_input_wm_arr = np.array(bin_input_wm)
print(bin_input_wm_arr)


######ks test input binned case 2. median
binned_array_med, b ,c = plt.hist(N_Sigma_i_med,bins=bins_med,alpha=0.4,histtype='stepfilled', normed = False, edgecolor = 'black', linewidth=0.8)
bin_no = len(binned_array_med)
print(len(binned_array_med))
print(bins_med)
central_value_med = []

for s in range(1, bin_no+1) :
	central_value_med.append(0.5*(bins_med[s-1] + bins_med[s]))
central_value_arr_med = np.array(central_value_med)
print(central_value_arr_med)

bin_input_med = []
for s in range(0, 37) :
	for i in range(0, int(binned_array_med[s])) :
		bin_input_med.append(central_value_arr_med[s])

bin_input_med_arr = np.array(bin_input_med)
#print(bin_input_med_arr)

#######################ks test for unbinned data 1. mean





print("ks test with mean as the central estimate for unbinned case...")

D_value_gaussian, p_value_gaussian = stats.kstest(N_Sigma_i_wm,'norm')
print(stats.kstest(N_Sigma_i_wm,'norm'))
D_value_cauchy, p_value_cauchy = stats.kstest(N_Sigma_i_wm,'cauchy')
print(stats.kstest(N_Sigma_i_wm,'cauchy'))
D_value_dexp, p_value_dexp = stats.kstest(N_Sigma_i_wm, 'laplace')
print(stats.kstest(N_Sigma_i_wm, 'laplace'))

#for q in range(1,100) :
D_value_t, p_value_t = stats.kstest(N_Sigma_i_wm,'t', args=(1, ))
print(str(stats.kstest(N_Sigma_i_wm,'t', args=(1,))))


print("for Gaussian, p value is " + str(p_value_gaussian))
print("for cauchy, p value is " + str(p_value_cauchy))
print("for double exponential, p value is " + str(p_value_dexp))
print("for student's t, p value is " + str(p_value_t))

#######################ks test for unbinned data 2. median

print("ks test with median as the central estimate for unbinned case...")
D_value_gaussian, p_value_gaussian = stats.kstest(N_Sigma_i_med,'norm')
print(stats.kstest(N_Sigma_i_med,'norm'))
D_value_cauchy, p_value_cauchy = stats.kstest(N_Sigma_i_med,'cauchy')
print(stats.kstest(N_Sigma_i_med,'cauchy'))
D_value_dexp, p_value_dexp = stats.kstest(N_Sigma_i_med, 'laplace')
print(stats.kstest(N_Sigma_i_med, 'laplace'))

#for q in range(1,100) :
D_value_t, p_value_t = stats.kstest(N_Sigma_i_med,'t', args=(2, ))
print(stats.kstest(N_Sigma_i_med,'t', args=(2, )))


print("for Gaussian, p value is " + str(p_value_gaussian))
print("for cauchy, p value is " + str(p_value_cauchy))
print("for double exponential, p value is " + str(p_value_dexp))
print("for student's t, p value is " + str(p_value_t))

#######################ks test for binned data 1. mean

print("ks test with mean as the central estimate for binned case...")
D_value_gaussian, p_value_gaussian = stats.kstest(bin_input_wm_arr,'norm')
print(stats.kstest(bin_input_wm,'norm'))
D_value_cauchy, p_value_cauchy = stats.kstest(bin_input_wm_arr,'cauchy')
print(stats.kstest(bin_input_wm,'cauchy'))
D_value_dexp, p_value_dexp = stats.kstest(bin_input_wm_arr, 'laplace')
print(stats.kstest(bin_input_wm, 'laplace'))



#for q in range(1,100) :
D_value_t, p_value_t = stats.kstest(bin_input_wm_arr,'t', args=(1, ))
print(stats.kstest(bin_input_wm,'t', args=(1, )))


print("for Gaussian, p value is " + str(p_value_gaussian))
print("for cauchy, p value is " + str(p_value_cauchy))
print("for double exponential, p value is " + str(p_value_dexp))
print("for student's t, p value is " + str(p_value_t))



#######################ks test for binned data 2. median

print("ks test with meadian as the central estimate for binned case...")
D_value_gaussian, p_value_gaussian = stats.kstest(bin_input_med_arr,'norm')
print(stats.kstest(bin_input_med,'norm'))
D_value_cauchy, p_value_cauchy = stats.kstest(bin_input_med_arr,'cauchy')
print(stats.kstest(bin_input_med,'cauchy'))
D_value_dexp, p_value_dexp = stats.kstest(bin_input_med_arr, 'laplace')
print(stats.kstest(bin_input_med, 'laplace'))

#for q in range(1,100) :
D_value_t, p_value_t = stats.kstest(bin_input_med_arr,'t', args=(1, 0, 1))
print(stats.kstest(bin_input_med,'t', args=(1, 0, 1)))


print("for Gaussian, p value is " + str(p_value_gaussian))
print("for cauchy, p value is " + str(p_value_cauchy))
print("for double exponential, p value is " + str(p_value_dexp))
print("for student's t, p value is " + str(p_value_t))





fig, ax = plt.subplots(nrows=2, ncols=2)
plt.tight_layout()

plt.subplot(2, 2, 1)
plt.hist(N_Sigma_i_wm,bins=bins_wm,alpha=0.4,histtype='stepfilled', normed = True, edgecolor = 'black', linewidth=0.8)
plt.xlabel("$N_σ$", fontsize = 15)
plt.ylabel("probability density", fontsize = 15)
plt.plot(pdf_x_wm,fitted_pdf_wm,"black", linestyle="dashed", linewidth=1.5)

plt.subplot(2, 2, 2)
plt.hist(N_sigma_symm_arr_wm,bins=np.arange(-10-0.5,10+0.5,0.5),alpha=0.4,histtype='stepfilled', normed = True, weights=0.5*np.ones(len(N_sigma_symm_arr_wm)), edgecolor = 'black', linewidth=0.8)
plt.xlabel("$|N_σ|$", fontsize = 15)
plt.ylabel("probability density", fontsize = 15)
plt.plot(pdf_x1_wm,fitted_pdf1_wm,"black",linestyle="dashed", linewidth=1.5)


plt.subplot(2, 2, 3)

plt.hist(N_Sigma_i_med,bins=bins_med,alpha=0.4,histtype='stepfilled', normed = True, edgecolor = 'black', linewidth=0.8)
plt.xlabel("$N_σ$", fontsize = 15)
plt.ylabel("probability density", fontsize = 15)
plt.plot(pdf_x_med,fitted_pdf_med,"black", linestyle="dashed", linewidth=1.5)

plt.subplot(2, 2, 4)
plt.hist(N_sigma_symm_arr_med,bins=np.arange(-10-0.5,10+0.5,0.5),alpha=0.4,histtype='stepfilled', normed = True, weights=0.5*np.ones(len(N_sigma_symm_arr_med)), edgecolor = 'black', linewidth=0.8)
plt.xlabel("$|N_σ|$", fontsize = 15)
plt.ylabel("probability density", fontsize = 15)
plt.plot(pdf_x1_med,fitted_pdf1_med,"black",linestyle="dashed", linewidth=1.5)

#plt.show()


















fig, ax = plt.subplots(nrows=3, ncols=2)
plt.tight_layout()


plt.subplot(3, 2, 1)
plt.hist(N_Sigma_i_wm,bins=bins_wm,alpha=0.4,histtype='stepfilled', normed = True, edgecolor = 'black', linewidth=0.8)
plt.xlabel("$N_σ$", fontsize = 14)
plt.ylabel("probability density", fontsize = 14)
parameters = stats.cauchy.fit(N_Sigma_i_wm)
pdf_x = np.linspace(np.amin(N_Sigma_i_wm) - 0.5,np.amax(N_Sigma_i_wm) + 0.5 ,500)
fitted_pdf_cauchy_wm = stats.cauchy.pdf(pdf_x,loc = parameters[0],scale = parameters[1])
plt.plot(pdf_x,fitted_pdf_cauchy_wm,"black", linestyle="dashed", linewidth=1.5)

plt.subplot(3, 2, 2)
plt.hist(N_Sigma_i_med,bins=bins_med,alpha=0.4,histtype='stepfilled', normed = True, edgecolor = 'black', linewidth=0.8)
plt.xlabel("$N_σ$", fontsize = 14)
#plt.ylabel("probability density")
parameters = stats.cauchy.fit(N_Sigma_i_med)
pdf_x = np.linspace(np.amin(N_Sigma_i_med) - 0.5,np.amax(N_Sigma_i_med) + 0.5 ,500)
fitted_pdf_cauchy_med = stats.cauchy.pdf(pdf_x,loc = parameters[0],scale = parameters[1])
plt.plot(pdf_x,fitted_pdf_cauchy_med,"black", linestyle="dashed", linewidth=1.5)

#plt.show()



plt.subplot(3, 2, 3)
plt.hist(N_Sigma_i_wm,bins=bins_wm,alpha=0.4,histtype='stepfilled', normed = True, edgecolor = 'black', linewidth=0.8)
plt.xlabel("$N_σ$", fontsize = 14)
plt.ylabel("probability density", fontsize = 14)
parameters = stats.laplace.fit(N_Sigma_i_wm)
pdf_x = np.linspace(np.amin(N_Sigma_i_wm) - 0.5,np.amax(N_Sigma_i_wm) + 0.5 ,500)
fitted_pdf_laplace_wm = stats.laplace.pdf(pdf_x,loc = parameters[0],scale = parameters[1])
plt.plot(pdf_x,fitted_pdf_laplace_wm,"black", linestyle="dashed", linewidth=1.5)

plt.subplot(3, 2, 4)
plt.hist(N_Sigma_i_med,bins=bins_med,alpha=0.4,histtype='stepfilled', normed = True, edgecolor = 'black', linewidth=0.8)
plt.xlabel("$N_σ$", fontsize = 14)
#plt.ylabel("probability density")
parameters = stats.laplace.fit(N_Sigma_i_med)
pdf_x = np.linspace(np.amin(N_Sigma_i_med) - 0.5,np.amax(N_Sigma_i_med) + 0.5 ,500)
fitted_pdf_laplace_med = stats.laplace.pdf(pdf_x,loc = parameters[0],scale = parameters[1])
plt.plot(pdf_x,fitted_pdf_laplace_med,"black", linestyle="dashed", linewidth=1.5)

#plt.show()


plt.subplot(3, 2, 5)
plt.hist(N_Sigma_i_wm,bins=bins_wm,alpha=0.4,histtype='stepfilled', normed = True, edgecolor = 'black', linewidth=0.8)
plt.xlabel("$N_σ$",fontsize = 14)
plt.ylabel("probability density", fontsize = 14)
parameters = stats.t.fit(N_Sigma_i_wm)
print(parameters)
pdf_x = np.linspace(np.amin(N_Sigma_i_wm) - 0.5,np.amax(N_Sigma_i_wm) + 0.5 ,500)
fitted_pdf_t_wm = stats.t.pdf(pdf_x,parameters[0],parameters[1], parameters[2])
plt.plot(pdf_x,fitted_pdf_t_wm,"black", linestyle="dashed", linewidth=1.5)

plt.subplot(3, 2, 6)
plt.hist(N_Sigma_i_med,bins=bins_med,alpha=0.4,histtype='stepfilled', normed = True, edgecolor = 'black', linewidth=0.8)
plt.xlabel("$N_σ$", fontsize = 14)
#plt.ylabel("probability density")
parameters = stats.t.fit(N_Sigma_i_med)
print(parameters)
pdf_x = np.linspace(np.amin(N_Sigma_i_med) - 0.5,np.amax(N_Sigma_i_med) + 0.5 ,500)
fitted_pdf_t_med = stats.t.pdf(pdf_x,parameters[0],parameters[1], parameters[2])
plt.plot(pdf_x,fitted_pdf_t_med,"black", linestyle="dashed", linewidth=1.5)

plt.show()
plt.savefig("image.png")

##################
fig, ax = plt.subplots(nrows=1, ncols=2)

plt.subplot(1, 2, 1)
plt.hist(N_Sigma_i_wm,bins=bins_wm,alpha=0.4,histtype='stepfilled', normed = False, edgecolor = 'black', linewidth=0.8)
plt.xlabel("$N_σ$")
plt.ylabel("probability density")
parameters = stats.weibull_min.fit(N_Sigma_i_wm)
pdf_x = np.linspace(np.amin(N_Sigma_i_wm) - 0.5,np.amax(N_Sigma_i_wm) + 0.5 ,500)
fitted_pdf = stats.weibull_min.pdf(pdf_x,parameters[0],parameters[1],parameters[2])
plt.plot(pdf_x,fitted_pdf,"black", linestyle="dashed", linewidth=1.5)

plt.subplot(1, 2, 2)
plt.hist(N_Sigma_i_med,bins=bins_med,alpha=0.4,histtype='stepfilled', normed = False, edgecolor = 'black', linewidth=0.8)
plt.xlabel("$N_σ$")
plt.ylabel("probability density")
parameters = stats.weibull_min.fit(N_Sigma_i_med)
pdf_x = np.linspace(np.amin(N_Sigma_i_med) - 0.5,np.amax(N_Sigma_i_med) + 0.5 ,500)
fitted_pdf = stats.weibull_min.pdf(pdf_x,parameters[0],parameters[1],parameters[2])
plt.plot(pdf_x,fitted_pdf,"black", linestyle="dashed", linewidth=1.5)

#plt.show()


'''

########################Chi square test
n_bins= int(np.sqrt(134))
print(n_bins)
N_sigma_chi_wm = N_sigma_symm_arr_wm[np.logical_and(N_sigma_symm_arr_wm<=10, N_sigma_symm_arr_wm>=-10)]
#print(N_sigma_chi_wm)
bin_size_chi_wm= (np.amax(N_sigma_chi_wm ) - np.amin(N_sigma_chi_wm))/n_bins

hist_wm, bin_edges_wm= np.histogram(N_sigma_chi_wm, bins=np.arange(np.amin(N_sigma_chi_wm ) ,np.amax(N_sigma_chi_wm), bin_size_chi_wm ) , normed = False)
print(hist_wm)
bin_centre_wm = []

for x in range(1,int(len(bin_edges_wm))) :
	bin_centre_wm.append(0.5*(bin_edges_wm[x-1] + bin_edges_wm[x]))

bin_center_arr_wm = np.array(bin_centre_wm)


N_sigma_chi_med = N_sigma_symm_arr_med[np.logical_and(N_sigma_symm_arr_med<=10, N_sigma_symm_arr_med>=-10)]

bin_size_chi_med= (np.amax(N_sigma_chi_med ) - np.amin(N_sigma_chi_med))/n_bins
hist_med, bin_edges_med= np.histogram(N_sigma_chi_med,bins=np.arange(np.amin(N_sigma_chi_wm ) ,np.amax(N_sigma_chi_wm) ,bin_size_chi_med), normed = False)
print(hist_med)
#print(bin_edges_med)

bin_centre_med = []

for x in range(1,int(len(bin_edges_med))) :
	bin_centre_med.append(0.5*(bin_edges_med[x-1] + bin_edges_med[x]))

bin_center_arr_med = np.array(bin_centre_med)


n_wm = 134*stats.norm.pdf(bin_center_arr_wm, 0, 1)
n_med = 134*stats.norm.pdf(bin_center_arr_med, 0, 1)
c_wm = 134*stats.cauchy.pdf(bin_center_arr_wm, 0, 1)
c_med = 134*stats.cauchy.pdf(bin_center_arr_med, 0, 1)
l_wm = 134*stats.laplace.pdf(bin_center_arr_wm, 0, 1)
l_med = 134*stats.laplace.pdf(bin_center_arr_med, 0, 1)
t_wm = 134*stats.t.pdf(bin_center_arr_wm,1, 0, 1)
t_med = 134*stats.t.pdf(bin_center_arr_med,1, 0, 1)

s=0
for d in range(1,10):
	s+=0.1
	k= 1/s
	hist_wm_scaled= hist_wm*k
	a= plt.scatter(bin_center_arr_wm, np.log10(hist_wm_scaled))
	b= plt.scatter(bin_center_arr_wm, np.log10(n_wm))
	plt.legend((a, b), ('term1', 'term2'))
	plt.show()

chi_sq_wm_norm = math.fsum((hist_wm - n_wm)*(hist_wm - n_wm)/n_wm)
chi_sq_med_norm = math.fsum((hist_med - n_med)*(hist_med - n_med)/n_med)

chi_sq_wm_cauchy = math.fsum((hist_wm - c_wm)*(hist_wm - c_wm)/c_wm)
chi_sq_med_cauchy = math.fsum((hist_med - c_med)*(hist_med - c_med)/c_med)

chi_sq_wm_laplace= math.fsum((hist_wm - l_wm)*(hist_wm - l_wm)/l_wm)
chi_sq_med_laplace = math.fsum((hist_med - l_med)*(hist_med - l_med)/l_med)

chi_sq_wm_t= math.fsum((hist_wm - t_wm)*(hist_wm - t_wm)/t_wm)
chi_sq_med_t = math.fsum((hist_med - t_med)*(hist_med - t_med)/t_med)


p_value_wm_norm = 1 - stats.chi2.cdf(x=chi_sq_wm_norm,  df=n_bins-1)
p_value_med_norm = 1 - stats.chi2.cdf(x=chi_sq_med_norm,  df=n_bins-1)
p_value_wm_cauchy = 1 - stats.chi2.cdf(x=chi_sq_wm_cauchy,  df=n_bins-1)
p_value_med_cauchy = 1 - stats.chi2.cdf(x=chi_sq_med_cauchy,  df=n_bins-1)
p_value_wm_laplace = 1 - stats.chi2.cdf(x=chi_sq_wm_laplace,  df=n_bins-1)

p_value_med_laplace = 1 - stats.chi2.cdf(x=chi_sq_med_laplace,  df=n_bins-1)
p_value_wm_t = 1 - stats.chi2.cdf(x=chi_sq_wm_t,  df=n_bins-2)
p_value_med_t = 1 - stats.chi2.cdf(x=chi_sq_med_t,  df=n_bins-2)

print(stats.chi2(n_bins -1).pdf(chi_sq_wm_norm))
print(stats.chi2(n_bins -1).pdf(chi_sq_med_norm))
print(stats.chi2(n_bins -1).pdf(chi_sq_wm_cauchy))
print(stats.chi2(n_bins -1).pdf(chi_sq_med_cauchy))
print(stats.chi2(n_bins -1).pdf(chi_sq_wm_laplace))
print(stats.chi2(n_bins -1).pdf(chi_sq_med_laplace))
print(stats.chi2(n_bins -2).pdf(chi_sq_wm_t))
print(stats.chi2(n_bins -2).pdf(chi_sq_med_t))


p_value_wm_norm = stats.chi2(n_bins-1).pdf(chi_sq_wm_norm)
p_value_med_norm =stats.chi2(n_bins-1).pdf(chi_sq_med_norm)
p_value_wm_cauchy =stats.chi2(n_bins-1).pdf(chi_sq_wm_cauchy)
p_value_med_cauchy =stats.chi2(n_bins-1)pdf(chi_sq_med_cauchy)


print("chi square test for mean and median values for each pdf...")
print(p_value_med_norm) 
print(p_value_wm_cauchy) 
print(p_value_med_cauchy) 
print(p_value_wm_laplace) 
print(p_value_med_laplace) 
print(p_value_wm_t) 
print(p_value_med_t) 
'''
