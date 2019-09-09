import pandas as pd
from matplotlib import pyplot as plt
import numpy as np

from relax import relaxation_fit, single_step_relaxation, two_step_relaxation, second_step_linear, linear, michaelis_menten

### CONSTANTS
color_set = [(0.2,0.6,0.2), (1,0.55,0.15)]
CONCENTRATION = 0.025
TIME = 18*60*60
CYCLE_NUMBER_INDEX = 56
# INITIAL_GUESS = [2000, 0.0001, 0.1, 100]
# INITIAL_GUESS = [2000, 0.0001, 100]
INITIAL_GUESS = 0.2, 100
### LOAD IN DATA

FILENAME = "/Users/benjaminbarad/Desktop/xlsx/FECN_Kinetic_20190528_154833.xlsx"

data = pd.read_excel(FILENAME, index_col=0, nrows=96, header=CYCLE_NUMBER_INDEX, skiprows=[CYCLE_NUMBER_INDEX+1]).transpose()
data = data.loc[data.index > 2000]  
### SORT DATA BY FUNCTION
cat = ["B","E"]
fl = ["F"]
no_enzyme = ["A"]
# wt = ["B","C"]
# fl = ["F","D","E"]
# no_enzyme = ["A", "G"]
std = [["H{}".format(i)] for i in range(1,7)]
conc = [2*(1/3)**i for i in range(11)] + [0]
# std_conc = [500 / 2**i for i in range(5)] + [0]
std_conc = [250, 100, 50, 20, 10, 0]

cat_results = []
fl_results = []
fig,ax = plt.subplots()
initial_rates_cat = []
initial_std_cat = []
initial_rates_fl = []
initial_std_fl = []

std_vals = []
std_stds = []
for index, concentration in enumerate(std_conc):
	standard_letters = std[index]
	standard = data[standard_letters].mean(axis=1)
	standard_std_set = data[standard_letters].std(axis=1)
	# values, covariances, y_calc = relaxation_fit(standard.index.values, standard, relaxation_function = linear, initial_guess=(-0.0002, 1.), sigma=standard_std_set)
	value = max(standard)-min(standard)
	standard_std = standard.std()
	# standard = standard.mean()
	print (standard)
	# standard_std = data[standard_letters].mean(axis=0).std(axis=2)
	std_vals.append(value)
	std_stds.append(standard.std())

adjuster, covariances, y_calc = relaxation_fit( std_vals, std_conc, relaxation_function=linear, initial_guess=(100, 100)) #, sigma=std_stds[1:11])
slope = adjuster[0]
intercept = adjuster[1]

fig1, ax1 = plt.subplots()
ax1.plot(std_vals, std_conc ,".")
ax1.plot(std_vals, y_calc )
fig1.savefig("standard_series.png")

for index, concentration in enumerate(conc):
	

	cat_letters = ["{0}{1}".format(i, index+1) for i in cat]
	scores = []
	for letter in cat_letters:
		trace = data[letter]
		difference = max(trace)-min(trace)
		scores.append(difference*slope+intercept)

	# cat_adjusted = np.mean(scores)
	
	# rates, covariances, y_calc = relaxation_fit(cat_adjusted.index.values, cat_adjusted, relaxation_function=linear, initial_guess = INITIAL_GUESS, sigma=cat_std, maxfev=30000)
	initial_rates_cat.append(np.mean(scores))
	initial_std_cat.append(np.std(scores))
	# ax.plot(cat_adjusted.index.values, cat_adjusted)
	# ax.plot(cat_adjusted.index.values, y_calc)

	fl_letters = ["{0}{1}".format(i, index+1) for i in fl]
	# fl_unadjusted = data[fl_letters].mean(axis=1)
	# fl_adjusted = fl_unadjusted - no_enzyme_control
	# fl_std = data[fl_letters].mean(axis=1)
	scores = []
	for letter in fl_letters:
		trace = data[letter]
		difference = max(trace)-min(trace)
		scores.append(difference*slope+intercept)
	# rates, covariances, y_calc = relaxation_fit(fl_adjusted.index.values, fl_adjusted, relaxation_function=linear, initial_guess = INITIAL_GUESS, sigma=fl_std, maxfev=30000)
	initial_rates_fl.append(np.mean(scores))
	print(scores)
	initial_std_fl.append(np.std(scores))
	# ax.plot(fl_adjusted.index.values, fl_adjusted)
	# ax.plot(fl_adjusted.index.values, y_calc)



# fig.savefig("adjusted.png")



### Check that the data load actually worked
fig2,ax2 = plt.subplots()
ax2.plot(data.loc[0:,"A1":"A12"])
fig2.savefig("pandas_sanity_check.png")


### Michaelis Menten Plot
fig3,ax3 = plt.subplots()

values, covar, y_calc = relaxation_fit(conc[1:11], initial_rates_cat[1:11], relaxation_function=michaelis_menten, initial_guess=(0.05,20), maxfev=30000)
print (values)
ax3.plot(conc[1:11], initial_rates_cat[1:11], '.', label="CatD", color=color_set[0]) #yerr=initial_std_wt[3:11]
ax3.plot(conc[1:11], y_calc, color=color_set[0])

values_fl, covar_fl, y_calc = relaxation_fit(conc[1:11], initial_rates_fl[1:11], relaxation_function=michaelis_menten, initial_guess=(0.05,20), maxfev=30000)
print(values_fl)
ax3.plot(conc[1:11], initial_rates_fl[1:11], '.', label = "FL", color=color_set[1]) # , yerr=initial_std_mut[3:11]
ax3.plot(conc[1:11], y_calc, color=color_set[1])

ax3.legend(loc=3)
fig3.savefig("MM_fecn.png")



print("CatD")
print("Vmax: {} ± {}".format(values[0]/TIME, np.sqrt(covar[0][0])/TIME))
print("kcat: {} ± {}".format(values[0]/CONCENTRATION/TIME, np.sqrt(covar[0][0])/CONCENTRATION/TIME))
print("Km: {} ± {}".format(values[1], np.sqrt(covar[1][1])))
cat_rel_act = values[0]/values[1]/CONCENTRATION/TIME
cat_rel_act_dev = cat_rel_act*np.sqrt(covar[0][0]/(values[0]**2)+covar[1][1]/(values[1]**2)-2*covar[0][1]/(values[0]*values[1]))
print(cat_rel_act, cat_rel_act_dev)

print("FL")
print("Vmax: {} ± {}".format(values_fl[0]/TIME, np.sqrt(covar_fl[0][0])/TIME))
print("kcat: {} ± {}".format(values_fl[0]/CONCENTRATION/TIME, np.sqrt(covar_fl[0][0])/CONCENTRATION/TIME))
print("Km: {} ± {}".format(values_fl[1], np.sqrt(covar_fl[1][1])))
fl_rel_act = values_fl[0]/values_fl[1]/CONCENTRATION/TIME
fl_rel_act_dev = fl_rel_act*np.sqrt(covar_fl[0][0]/(values_fl[0]**2)+covar_fl[1][1]/(values_fl[1]**2)-2*covar_fl[0][1]/(values_fl[0]*values_fl[1]))
print(fl_rel_act, fl_rel_act_dev)

## Bar Plot
fig4,ax4=plt.subplots()
x = [0,1]
y = [cat_rel_act, fl_rel_act]
y_tick_labels = ["CatD", "FL"]
# y = [i/wt_rel_act for i in y]
yerr = [cat_rel_act_dev, fl_rel_act_dev]
# yerr = [i/wt_rel_act for i in yerr]
color_set = [(0.2,0.6,0.2), (1,0.55,0.15)]
ax4.set_xticks(x)
ax4.set_xticklabels(y_tick_labels)
ax4.bar(x, y, 0.8, yerr=yerr, color=color_set)
ax4.set_ylabel(r"kcat/Km $(s•[S])^{-1}$")
fig4.savefig("1C.png")

kcat_km_ratio = y[1]/y[0]
kcat_km_ratio_error = kcat_km_ratio * np.sqrt((yerr[0]/y[0])**2 + (yerr[1]/y[1])**2)

fig6,ax6=plt.subplots()
x = [0,1]
y = [values[1], values_fl[1]]
y_tick_labels = ["CatD", "FL"]
# y = [i/wt_rel_act for i in y]
yerr = [np.sqrt(covar[1][1]), np.sqrt(covar_fl[1][1])]
# yerr = [i/wt_rel_act for i in yerr]
color_set = [(0.2,0.6,0.2), (1,0.55,0.15)]
ax6.set_xticks(x)
ax6.set_xticklabels(y_tick_labels)
ax6.bar(x, y, 0.8, yerr=yerr, color=color_set)
ax6.set_ylabel(r"Km $([S])^{-1}$")
plt.tight_layout()
fig6.savefig("1C_Km.png")

km_ratio = y[1]/y[0]
km_ratio_error = kcat_km_ratio * np.sqrt((yerr[0]/y[0])**2 + (yerr[1]/y[1])**2)

fig5,ax5=plt.subplots()
x = [0,1]
y = [values[0]/CONCENTRATION/TIME, values_fl[0]/CONCENTRATION/TIME]
y_tick_labels = ["CatD", "FL"]
# y = [i/wt_rel_act for i in y]
yerr = [np.sqrt(covar[0][0])/CONCENTRATION/TIME, np.sqrt(covar_fl[0][0])/CONCENTRATION/TIME]
# yerr = [i/wt_rel_act for i in yerr]
color_set = [(0.2,0.6,0.2), (1,0.55,0.15)]
ax5.set_xticks(x)
ax5.set_xticklabels(y_tick_labels)
ax5.bar(x, y, 0.8, yerr=yerr, color=color_set)
ax5.set_ylabel(r"kcat $(s)^{-1}$")
plt.tight_layout()
fig5.savefig("1C_kcat.png")


kcat_ratio = y[1]/y[0]
kcat_ratio_error = kcat_km_ratio * np.sqrt((yerr[0]/y[0])**2 + (yerr[1]/y[1])**2)


fig7, ax7 = plt.subplots()
y = [kcat_ratio, km_ratio, kcat_km_ratio]
yerr = [kcat_ratio_error, km_ratio_error, kcat_km_ratio_error]
x = [0,1,2]
x_tick_labels = ["kcat", "Km", "kcat/Km"]
ax7.set_ylabel("Ratio of Activity (FL/CatD)")
ax7.set_xlabel("Component of Activity")
color_set = [(0, 114./256, 178./256), (204./256, 121./256, 167./256), (230./256, 159./256, 0)]
ax7.set_xticks(x)
ax7.set_xticklabels(x_tick_labels)
ax7.bar(x,y,0.8, yerr=yerr, color=color_set)
plt.tight_layout()
fig7.savefig("1C_ratios.png")




