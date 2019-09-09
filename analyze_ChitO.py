import pandas as pd
from matplotlib import pyplot as plt
import numpy as np

from relax import relaxation_fit, single_step_relaxation, two_step_relaxation, second_step_linear, linear, michaelis_menten

### CONSTANTS
color_set = [(0.2,0.6,0.2), (1,0.55,0.15)]
CONCENTRATION = 0.001
CYCLE_NUMBER_INDEX = 65 #quantared
# CYCLE_NUMBER_INDEX = 73 #quantablu
# INITIAL_GUESS = [2000, 0.0001, 0.1, 100]
INITIAL_GUESS = [2000, 0.0001, 100]
# INITIAL_GUESS = [0.2, 100]
### LOAD IN DATA

FILENAME = "/Users/benjaminbarad/Desktop/xlsx/ChitO_Test_1_20190425_092146.xlsx"
# FILENAME = "/Users/benjaminbarad/Desktop/xlsx/ChitO_Test_1_20190521_081704.xlsx"
# FILENAME = "/Users/benjaminbarad/Desktop/xlsx/ChitO_QuantaBlu_20190321_051503.xlsx"

data = pd.read_excel(FILENAME, index_col=0, nrows=96, header=CYCLE_NUMBER_INDEX, skiprows=[CYCLE_NUMBER_INDEX+1]).transpose()
# data = data.loc[data.index > 500]  
data = data.loc[data.index<30000]
### SORT DATA BY FUNCTION
wt = ["B","C"]
mut = ["E","D"]
no_enzyme = ["G"]
# wt = ["B","C"]
# mut = ["F","D","E"]
# no_enzyme = ["A", "G"]
std = [["H{}".format(i+6)] for i in range(1,7)] #"H{}".format(i), 
conc = [0.25*(2/3)**i for i in range(11)] + [0]
std_conc = [50 / 2**i for i in range(5)] + [0]


wt_results = []
mut_results = []
fig,ax = plt.subplots()
initial_rates_wt = []
initial_std_wt = []
initial_rates_mut = []
initial_std_mut = []

std_vals = []
std_stds = []
for index, concentration in enumerate(std_conc):
	standard_letters = std[index]
	standard = data[standard_letters].loc[data.index>10000].loc[data.index<15000].mean(axis=1)
	standard_std = standard.std()
	standard = standard.mean()
	# standard_std = data[standard_letters].mean(axis=0).std(axis=2)
	std_vals.append(standard)
	std_stds.append(standard_std)

adjuster, covariances, y_calc = relaxation_fit(std_vals[1:], std_conc[1:], relaxation_function=linear, initial_guess=(100, 100), sigma=std_stds[1:])
# adjuster = [ 0.00729535, -1.2748655 ] # from another run, since the standards failed on this one
slope = adjuster[0]
intercept = adjuster[1]
print(adjuster)


fig1, ax1 = plt.subplots()
ax1.plot(std_conc, std_vals ,".")
ax1.plot(y_calc, std_vals[1:])
ax1.set_xlabel(r"Standard Concentration ($\mu M$)")
ax1.set_ylabel("Fluorescence (RFU)")
plt.tight_layout()
fig1.savefig("standard_series.png")


for index, concentration in enumerate(conc):
	

	no_enzyme_letters = ["{0}{1}".format(i, index+1) for i in no_enzyme]
	no_enzyme_control = data[no_enzyme_letters].mean(axis=1)
	# ax.plot(no_enzyme_control.index.values, no_enzyme_control)
	
	wt_letters = ["{0}{1}".format(i, index+1) for i in wt]
	wt_unadjusted = data[wt_letters].mean(axis=1)
	wt_adjusted = wt_unadjusted - no_enzyme_control
	wt_std = data[wt_letters].mean(axis=1)
	
	rates, covariances, y_calc = relaxation_fit(wt_adjusted.index.values, wt_adjusted, relaxation_function=single_step_relaxation, initial_guess = INITIAL_GUESS, maxfev=30000)
	initial_rates_wt.append(rates[2]*slope)
	initial_std_wt.append(np.sqrt(covariances[0][0]/(rates[0]**2) + covariances[1][1]/(rates[1]**2) + 2*covariances[0][1]/(rates[0]*rates[1]))*rates[0]*rates[1]*slope)

	mut_letters = ["{0}{1}".format(i, index+1) for i in mut]
	mut_unadjusted = data[mut_letters].mean(axis=1)
	mut_adjusted = mut_unadjusted - no_enzyme_control
	mut_std = data[mut_letters].mean(axis=1)
	
	rates, covariances, y_calc = relaxation_fit(mut_adjusted.index.values, mut_adjusted, relaxation_function=single_step_relaxation, initial_guess = INITIAL_GUESS, maxfev=30000)
	initial_rates_mut.append(rates[2]*slope) #*rates[1]
	print(rates[0])
	initial_std_mut.append(np.sqrt(covariances[0][0]/(rates[0]**2) + covariances[1][1]/(rates[1]**2) + 2*covariances[0][1]/(rates[0]*rates[1]))*rates[0]*rates[1]*slope)
	if index>5:
		ax.plot(wt_adjusted.index.values, wt_adjusted)
		ax.plot(mut_adjusted.index.values, mut_adjusted)
		ax.plot(mut_adjusted.index.values, y_calc)

	# ax.plot(mut_adjusted.index.values, mut_adjusted)
	# ax.plot(mut_adjusted.index.values, y_calc)


ax.set_ylabel("Fluorescence (RFU)")
ax.set_xlabel("Time (s)")
plt.tight_layout()
fig.savefig("adjusted.png")



### Check that the data load actually worked
fig2,ax2 = plt.subplots()
ax2.plot(data.loc[0:,"H1":"H12"])
fig2.savefig("pandas_sanity_check.png")

### Michaelis Menten Plot
fig3,ax3 = plt.subplots()

values, covar, y_calc = relaxation_fit(conc[4:11], initial_rates_wt[4:11], relaxation_function=michaelis_menten, initial_guess=(0.05,20), maxfev=30000) #, absolute_sigma=True)
print (values)
ax3.plot(conc[4:11], initial_rates_wt[4:11], '.', label="WT", color=color_set[0]) #yerr=initial_std_wt[3:11]
ax3.plot(conc[4:11], y_calc, color=color_set[0])

values_mut, covar_mut, y_calc = relaxation_fit(conc[4:11], initial_rates_mut[4:11],  relaxation_function=michaelis_menten, initial_guess=(0.05,20), maxfev=30000) #, absolute_sigma=True)
print(values_mut)
ax3.plot(conc[4:11], initial_rates_mut[4:11], '.', label = "Mut", color=color_set[1]) # , yerr=initial_std_mut[3:11]
ax3.plot(conc[4:11], y_calc, color=color_set[1])

ax3.legend(loc=4)
ax3.set_ylabel(r"Calculated Initial Rate ($\mu M/s$)")
ax3.set_xlabel("Chitin Concentration (% w/v)")
plt.tight_layout()
fig3.savefig("MM_chito.png")



print("WT")
print("Vmax: {} ± {}".format(values[0], np.sqrt(covar[0][0])))
print("kcat: {} ± {}".format(values[0]/CONCENTRATION, np.sqrt(covar[0][0])/CONCENTRATION))
print("Km: {} ± {}".format(values[1], np.sqrt(covar[1][1])))
wt_rel_act = values[0]/values[1]/CONCENTRATION
wt_rel_act_dev = wt_rel_act*np.sqrt(covar[0][0]/(values[0]**2)+covar[1][1]/(values[1]**2)-2*covar[0][1]/(values[0]*values[1]))
print(wt_rel_act, wt_rel_act_dev)

print("Mut")
print("Vmax: {} ± {}".format(values_mut[0], np.sqrt(covar_mut[0][0])))
print("kcat: {} ± {}".format(values_mut[0]/CONCENTRATION, np.sqrt(covar_mut[0][0])/CONCENTRATION))
print("Km: {} ± {}".format(values_mut[1], np.sqrt(covar_mut[1][1])))
mut_rel_act = values_mut[0]/values_mut[1]/CONCENTRATION
mut_rel_act_dev = mut_rel_act*np.sqrt(covar_mut[0][0]/(values_mut[0]**2)+covar_mut[1][1]/(values_mut[1]**2)-2*covar_mut[0][1]/(values_mut[0]*values_mut[1]))
print(mut_rel_act, mut_rel_act_dev)

## Bar Plot
fig4,ax4=plt.subplots()
x = [0,1]
y = [wt_rel_act, mut_rel_act]
y_tick_labels = ["AMCase", "Mut"]
# y = [i/wt_rel_act for i in y]
yerr = [wt_rel_act_dev, mut_rel_act_dev]
# yerr = [i/wt_rel_act for i in yerr]
color_set = [(0.2,0.6,0.2), (1,0.55,0.15)]
ax4.set_xticks(x)
ax4.set_xticklabels(y_tick_labels)
ax4.bar(x, y, 0.8, yerr=yerr, color=color_set)
ax4.set_ylabel(r"kcat/Km $(s•\mu M)^{-1}$")
fig4.savefig("ChitO_Ratio.png")
