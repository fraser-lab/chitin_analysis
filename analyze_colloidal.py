import pandas as pd
from matplotlib import pyplot as plt
import numpy as np

from relax import relaxation_fit, single_step_relaxation, two_step_relaxation, second_step_linear, linear, michaelis_menten

### CONSTANTS
color_set = [(0.2,0.6,0.2), (1,0.55,0.15)]
CONCENTRATION = 0.005
CYCLE_NUMBER_INDEX = 54
# INITIAL_GUESS = [2000, 0.0001, 0.1, 100]
INITIAL_GUESS = [-1, 0.00001, 1]
# INITIAL_GUESS = -0.2, 1
### LOAD IN DATA

FILENAME = "/Users/benjaminbarad/Desktop/xlsx/Colloidal_Clearance_20190430_120255.xlsx"

data = pd.read_excel(FILENAME, index_col=0, nrows=96, header=CYCLE_NUMBER_INDEX, skiprows=[CYCLE_NUMBER_INDEX+1]).transpose()
data = data.loc[data.index > 5000]  
# data = data.loc[data.index < 50000]
### SORT DATA BY FUNCTION
fl = ["E", "F"]
cat = ["B", "C"]
no_enzyme = ["A"]
# wt = ["B","C"]
# fl = ["F","D","E"]
# no_enzyme = ["A", "G"]
# std = [["H{}".format(i)] for i in range(1,7)]
conc = [2*(1/2)**i for i in range(11)] + [0]
# std_conc = [500 / 2**i for i in range(5)] + [0]
# std_conc = [250, 100, 50, 20, 10, 0]

cat_results = []
fl_results = []
fig,ax = plt.subplots()
initial_rates_cat = []
initial_std_cat = []
initial_rates_fl = []
initial_std_fl = []

std_vals = []
std_stds = []
for index, concentration in enumerate(conc):
	standard_set = data["A{}".format(index+1)]
	standard = standard_set.mean()
	standard_std = standard_set.std()
	print (standard, standard_std)
	# standard_std = data[standard_letters].mean(axis=0).std(axis=2)
	std_vals.append(standard)
	std_stds.append(standard_std)

adjuster, covariances, y_calc = relaxation_fit( std_vals[1:11], conc[1:11], sigma=std_stds[1:11], relaxation_function=linear, initial_guess=(100, 100)) #, sigma=std_stds[1:11])
slope = adjuster[0]
intercept = adjuster[1]

fig1, ax1 = plt.subplots()
ax1.plot(std_vals[:-1], conc[:-1] ,".")
ax1.plot(std_vals[1:11], y_calc )
fig1.savefig("standard_series.png")

for index, concentration in enumerate(conc):
	no_enzyme_letters = ["{0}{1}".format(i, index+1) for i in no_enzyme]
	no_enzyme_control = data[no_enzyme_letters].mean(axis=1)

	cat_letters = ["{0}{1}".format(i, index+1) for i in cat]
	cat_unadjusted = data[cat_letters].mean(axis=1)
	# cat_adjusted =  no_enzyme_control - cat_unadjusted 
	cat_adjusted = cat_unadjusted.iloc[0] - cat_unadjusted
	print (cat_adjusted)
	cat_std = data[cat_letters].std(axis=1)

	rates, covariances, y_calc = relaxation_fit(cat_adjusted.index.values, cat_adjusted, relaxation_function=single_step_relaxation, initial_guess = INITIAL_GUESS, maxfev=30000)
	initial_rates_cat.append(rates[1]*rates[0])
	initial_std_cat.append(np.sqrt(covariances[0][0]/rates[0]**2+covariances[1][1]/rates[1]**2 + 2*covariances[1][0]/(rates[1]*rates[0]))*rates[0]*rates[1])


	fl_letters = ["{0}{1}".format(i, index+1) for i in fl]
	fl_unadjusted = data[fl_letters].mean(axis=1)
	fl_adjusted = fl_unadjusted.iloc[0] - fl_unadjusted

	rates, covariances, y_calc = relaxation_fit(fl_adjusted.index.values, fl_adjusted, relaxation_function=single_step_relaxation, initial_guess = INITIAL_GUESS, maxfev=30000)
	initial_rates_fl.append(rates[1]*rates[0])
	initial_std_fl.append(np.sqrt(covariances[0][0]/rates[0]**2+covariances[1][1]/rates[1]**2 + 2*covariances[1][0]/(rates[1]*rates[0]))*rates[0]*rates[1])

	if index > 1:
		ax.plot(fl_adjusted.index.values, fl_adjusted,'.')
		ax.plot(fl_adjusted.index.values, y_calc)

fig.savefig("adjusted.png")



### Check that the data load actually worked
fig2,ax2 = plt.subplots()
ax2.plot(data.loc[0:,"A1":"A12"])
fig2.savefig("pandas_sanity_check.png")


### Michaelis Menten Plot
fig3,ax3 = plt.subplots()

values, covar, y_calc = relaxation_fit(conc[2:9], initial_rates_cat[2:9], sigma=initial_std_fl[2:9], relaxation_function=michaelis_menten, initial_guess=(0.05,20), maxfev=30000)
print (values)
ax3.plot(conc[2:9], initial_rates_cat[2:9], '.', label="CatD", color=color_set[0]) #yerr=initial_std_wt[3:11]
ax3.plot(conc[2:9], y_calc, color=color_set[0])

values_fl, covar_fl, y_calc = relaxation_fit(conc[2:9], initial_rates_fl[2:9], sigma=initial_std_fl[2:9], relaxation_function=michaelis_menten, initial_guess=(0.05,20), maxfev=30000)
print(values_fl)
ax3.plot(conc[2:9], initial_rates_fl[2:9], '.', label = "FL", color=color_set[1]) # , yerr=initial_std_mut[3:11]
ax3.plot(conc[2:9], y_calc, color=color_set[1])

ax3.legend(loc=3)
fig3.savefig("MM_colloid.png")



print("CatD")
print("Vmax: {} ± {}".format(values[0]*slope, np.sqrt(covar[0][0])))
print("kcat: {} ± {}".format(values[0]/CONCENTRATION*slope, np.sqrt(covar[0][0])/CONCENTRATION))
print("Km: {} ± {}".format(values[1], np.sqrt(covar[1][1])))
cat_rel_act = values[0]/values[1]/CONCENTRATION
cat_rel_act_dev = cat_rel_act*np.sqrt(covar[0][0]/(values[0]**2)+covar[1][1]/(values[1]**2)+2*covar[0][1]/(values[0]*values[1]))
print(cat_rel_act, cat_rel_act_dev)

print("FL")
print("Vmax: {} ± {}".format(values_fl[0]*slope, np.sqrt(covar_fl[0][0])*slope))
print("kcat: {} ± {}".format(values_fl[0]/CONCENTRATION*slope, np.sqrt(covar_fl[0][0])/CONCENTRATION))
print("Km: {} ± {}".format(values_fl[1], np.sqrt(covar_fl[1][1])))
fl_rel_act = values_fl[0]/values_fl[1]/CONCENTRATION
fl_rel_act_dev = fl_rel_act*np.sqrt(covar_fl[0][0]/(values_fl[0]**2)+covar_fl[1][1]/(values_fl[1]**2)+2*covar_fl[0][1]/(values_fl[0]*values_fl[1]))
print(fl_rel_act, fl_rel_act_dev)

## Bar Plot
fig4,ax4=plt.subplots()
x = [0,1]
y = [cat_rel_act*slope, fl_rel_act*slope]
y_tick_labels = ["CatD", "FL"]
# y = [i/wt_rel_act for i in y]
yerr = [cat_rel_act_dev*slope, fl_rel_act_dev*slope]
# yerr = [i/wt_rel_act for i in yerr]
color_set = [(0.2,0.6,0.2), (1,0.55,0.15)]
ax4.set_xticks(x)
ax4.set_xticklabels(y_tick_labels)
ax4.bar(x, y, 0.8, yerr=yerr, color=color_set)
ax4.set_ylabel(r"kcat/Km $(s•\mu M)^{-1}$")
plt.tight_layout()
fig4.savefig("1B.png")
