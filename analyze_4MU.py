import pandas as pd
from matplotlib import pyplot as plt
import numpy as np

from relax import relaxation_fit, single_step_relaxation, two_step_relaxation, second_step_linear, linear, michaelis_menten

### CONSTANTS
CONCENTRATION = 0.05 # in uM - final enzyme concentration
CYCLE_NUMBER_INDEX = 59
# INITIAL_GUESS = [2000, 0.0001, 0.1, 100]
INITIAL_GUESS = [2000, 0.0001, 1000]
# INITIAL_GUESS = 0.2, 100
### LOAD IN DATA

# FILENAME = "/Users/benjaminbarad/Desktop/xlsx/4MU_Kinetic_Purified_20190410_205440.xlsx"
# FILENAME = "/Users/benjaminbarad/Desktop/xlsx/4MU_Kinetic_Purified (Modified)_20190403_125553_CatD_vs_N.xlsx"
FILENAME = "/Users/benjaminbarad/Desktop/xlsx/4MU_Kinetic_Purified_20190404_190520.xlsx"

data = pd.read_excel(FILENAME, index_col=0, nrows=96, header=CYCLE_NUMBER_INDEX, skiprows=[CYCLE_NUMBER_INDEX+1]).transpose()
# data = data.loc[data.index > 500]
data = data.loc[data.index<2000]
### SORT DATA BY FUNCTION
wt = ["B","A", "C"]
mut1 = ["E", "D", "F"]
mut2 = ["F", "E"]
std = ["H"]
conc = [108.5*(2/3)**i for i in range(11)] + [0]
std_conc = [177.4 *(1/2)**i for i in range(11)] + [0]

std_vals = []
std_stds = []
for index, concentration in enumerate(std_conc):
	standard_letters = ["{0}{1}".format(i, index+1) for i in std]
	standard = data[standard_letters].mean(axis=1)
	standard_std = standard.std()
	standard = standard.mean()
	print (standard)
	# standard_std = data[standard_letters].mean(axis=0).std(axis=2)
	std_vals.append(standard)
	std_stds.append(standard_std)
adjuster, covariances, y_calc = relaxation_fit(std_vals[2:11], std_conc[2:11], relaxation_function=linear, initial_guess=(100, 100), sigma=std_stds[2:11])
slope = adjuster[0]
intercept = adjuster[1]

fig1, ax1 = plt.subplots()
ax1.plot(std_vals[2:11], std_conc[2:11] ,".")
ax1.plot(std_vals[2:11], y_calc )
fig1.savefig("standard_series.png")

wt_results = []
mut_results = []
fig,ax = plt.subplots()
initial_rates_wt = []
initial_std_wt = []
initial_rates_mut = []
initial_std_mut = []
initial_rates_mut_2 =[]
initial_std_mut_2 = []
for index, concentration in enumerate(conc):


	# no_enzyme_letters = ["{0}{1}".format(i, index+1) for i in no_enzyme]
	# no_enzyme_control = data[no_enzyme_letters].mean(axis=1)

	wt_letters = ["{0}{1}".format(i, index+1) for i in wt]
	wt_unadjusted = data[wt_letters].mean(axis=1)
	wt_adjusted = wt_unadjusted
	print (wt_adjusted)
	wt_std = data[wt_letters].std(axis=1)

	rates, covariances, y_calc = relaxation_fit(wt_adjusted.index.values, wt_adjusted, relaxation_function=single_step_relaxation, initial_guess = INITIAL_GUESS, maxfev=30000) #, sigma=wt_std, absolute_sigma=True
	initial_rates_wt.append(rates[0]*rates[1]*slope)
	initial_std_wt.append(np.sqrt(covariances[0][0]/(rates[0]**2)+ covariances[1][1]/(rates[1]**2) + 2*covariances[0][1]/(rates[0]*rates[1]))*rates[0]*rates[1]*slope)
	# ax.plot(wt_adjusted.index.values, wt_adjusted)
	# ax.plot(wt_adjusted.index.values, y_calc)

	mut_letters = ["{0}{1}".format(i, index+1) for i in mut1]
	mut_unadjusted = data[mut_letters].mean(axis=1)
	mut_adjusted = mut_unadjusted
	mut_std = data[mut_letters].std(axis=1)
	rates, covariances, y_calc = relaxation_fit(mut_adjusted.index.values, mut_adjusted,relaxation_function=single_step_relaxation, initial_guess = INITIAL_GUESS, maxfev=30000)
	initial_rates_mut.append(rates[0]*rates[1]*slope)
	print(rates[0]*rates[1]*slope)
	initial_std_mut.append(np.sqrt(covariances[0][0]/(rates[0]**2)+ covariances[1][1]/(rates[1]**2) + 2*covariances[0][1]/(rates[0]*rates[1]))*rates[0]*rates[1]*slope)

	if index > 2:
		ax.plot(mut_adjusted.index.values, mut_adjusted,'.')
		ax.plot(mut_adjusted.index.values, y_calc)

	mut_2_letters = ["{0}{1}".format(i, index+1) for i in mut2]
	mut_2_unadjusted = data[mut_2_letters].mean(axis=1)
	mut_2_adjusted = mut_2_unadjusted
	mut_2_std = data[mut_2_letters].std(axis=1)
	rates, covariances, y_calc = relaxation_fit(mut_2_adjusted.index.values, mut_2_adjusted, sigma=mut_2_std, absolute_sigma=True, relaxation_function=single_step_relaxation, initial_guess = INITIAL_GUESS, maxfev=30000)
	initial_rates_mut_2.append(rates[0]*rates[1]*slope)
	initial_std_mut_2.append((np.sqrt(covariances[0][0]/(rates[0]**2)+ covariances[1][1]/(rates[1]**2) + 2*covariances[0][1]/(rates[0]*rates[1]))*rates[0]*rates[1]*slope))







# fig.savefig("adjusted.png")


### Check that the data load actually worked
fig2,ax2 = plt.subplots()
ax2.plot(data.loc[0:,"A1":"A12"])
fig2.savefig("pandas_sanity_check.png")

fig3,ax3 = plt.subplots()
values, covar, y_calc = relaxation_fit(conc[:9], initial_rates_wt[:9], relaxation_function=michaelis_menten, initial_guess=(0.05,20), maxfev=30000)
print (values)
ax3.plot(conc[:9], initial_rates_wt[:9], '.', label="WT") #yerr=initial_std_wt[3:11]
ax3.plot(conc[:9], y_calc)

values_mut_1, covar_mut_1, y_calc = relaxation_fit(conc[:9], initial_rates_mut[:9],relaxation_function=michaelis_menten, initial_guess=(0.5,20), maxfev=30000)
print(values_mut_1)
ax3.plot(conc[:9], initial_rates_mut[:9], '.', label = "Mutant 1") # , yerr=initial_std_mut[3:11]
ax3.plot(conc[:9], y_calc)

values_mut_2, covar_mut_2, y_calc = relaxation_fit(conc[1:11], initial_rates_mut_2[1:11], sigma=initial_std_mut_2[1:11],relaxation_function=michaelis_menten, initial_guess=(0.5,20), maxfev=30000)
# print(values_mut_2)
# ax3.plot(conc[1:11], initial_rates_mut_2[1:11], '.', label="Mutant 2")
# ax3.plot(conc[1:11], y_calc)

ax3.legend(loc=3)
fig3.savefig("MM.png")

print("Wild Type")
print("Vmax: {} ± {}".format(values[0], np.sqrt(covar[0][0])))
print("kcat: {} ± {}".format(values[0]/CONCENTRATION, np.sqrt(covar[0][0])/CONCENTRATION))
print("Km: {} ± {}".format(values[1], np.sqrt(covar[1][1])))
wt_rel_act = values[0]/values[1]/CONCENTRATION
wt_rel_act_dev = wt_rel_act*np.sqrt(covar[0][0]/(values[0]**2)+covar[1][1]/(values[1]**2)-2*covar[0][1]/(values[0]*values[1]))
print(wt_rel_act, wt_rel_act_dev)

print("Mutant A")
print("Vmax: {} ± {}".format(values_mut_1[0], np.sqrt(covar_mut_1[0][0])))
print("kcat: {} ± {}".format(values_mut_1[0]/CONCENTRATION, np.sqrt(covar_mut_1[0][0])/CONCENTRATION))
print("Km: {} ± {}".format(values_mut_1[1], np.sqrt(covar_mut_1[1][1])))
mut_rel_act = values_mut_1[0]/values_mut_1[1]/CONCENTRATION
mut_rel_act_dev = mut_rel_act*np.sqrt(covar_mut_1[0][0]/(values_mut_1[0]**2)+covar_mut_1[1][1]/(values_mut_1[1]**2)-2*covar_mut_1[0][1]/(values_mut_1[0]*values_mut_1[1]))
print(mut_rel_act, mut_rel_act_dev)


print("Mutant B")
print("Vmax: {} ± {}".format(values_mut_2[0], np.sqrt(covar_mut_2[0][0])))
print("kcat: {} ± {}".format(values_mut_2[0]/CONCENTRATION, np.sqrt(covar_mut_2[0][0])/CONCENTRATION))
print("Km: {} ± {}".format(values_mut_2[1], np.sqrt(covar_mut_2[1][1])))
mut_2_rel_act = values_mut_2[0]/values_mut_2[1]/CONCENTRATION
mut_2_rel_act_dev = mut_2_rel_act*np.sqrt(covar_mut_2[0][0]/(values_mut_2[0]**2)+covar_mut_2[1][1]/(values_mut_2[1]**2)-2*covar_mut_2[0][1]/(values_mut_2[0]*values_mut_2[1]))
print(mut_2_rel_act, mut_2_rel_act_dev)


fig4,ax4=plt.subplots()
x = [0,1,2]
y = [wt_rel_act, mut_rel_act, mut_2_rel_act]
y_tick_labels = ["WT", "A239T/L364Q", "V246A"]
# y = [i/wt_rel_act for i in y]
yerr = [wt_rel_act_dev, mut_rel_act_dev, mut_2_rel_act_dev]
# yerr = [i/wt_rel_act for i in yerr]
color_set = [(0.2,0.6,0.2), (1,0.55,0.15), (0.5,1,1)]
ax4.set_xticks(x)
ax4.set_xticklabels(y_tick_labels)
ax4.bar(x, y, 0.8, yerr=yerr, color=color_set)
ax4.set_ylabel(r"kcat/Km $(s•[S])^{-1}$")
fig4.savefig("3A.png")
