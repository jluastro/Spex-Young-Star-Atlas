# map Teff anf Spt Type
# using given values - interpolate rest
# Corey Mutnik 6/22/16

import scipy.optimize as optimized
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
from astropy.io import ascii
import numpy as np


# calculated EW values
spt_table = ascii.read('Teff.lumclassV.Astrophysical.Quantities.Ed4.pg151.dat')

teff = spt_table['T_eff[K]']
err_teff = spt_table['uncertainty[K]']

spt_type = list(spt_table['Spt_type'])
breakup_spt_type = list(map(list, spt_type))

mapped_val = []
for k in range(len(breakup_spt_type)):
	val = 0
	if breakup_spt_type[k][0] == 'B':
		val += 0
	elif breakup_spt_type[k][0] == 'A':
		val += 10
	elif breakup_spt_type[k][0] == 'F':
		val += 20
	elif breakup_spt_type[k][0] == 'G':
		val += 30
	elif breakup_spt_type[k][0] == 'K':
		val += 40
	elif breakup_spt_type[k][0] == 'M':
		val += 50
	else:
		val += 1000
		print 'what spt type is this:', breakup_spt_type[k][0],' index: ', k
	#print breakup_spt_type[k], val
	
	if breakup_spt_type[k][1] == '1':
		val += 1
	elif breakup_spt_type[k][1] == '2':
		val += 2
	elif breakup_spt_type[k][1] == '3':
		val += 3
	elif breakup_spt_type[k][1] == '4':
		val += 4
	elif breakup_spt_type[k][1] == '5':
		val += 5
	elif breakup_spt_type[k][1] == '6':
		val += 6
	elif breakup_spt_type[k][1] == '7':
		val += 7
	elif breakup_spt_type[k][1] == '8':
		val += 8
	elif breakup_spt_type[k][1] == '9':
		val += 9
	else:
		val += 0
	#print breakup_spt_type[k], val

	mapped_val.append(val)

'''
# function to fit
def func(x, a, b, c, d, e, f):
	return a + b*x + c*x**2 + d*x**3 + e*x**4 + f*x**5
# initial guesses
a,b,c,d,e,f = 1.,1.,1.,1.,1.,1.
param, covar = optimized.curve_fit(func, mapped_val, teff, p0=[a,b,c,d,e,f])
y_fit = func(x_fit, param[0], param[1], param[2], param[3], param[4], param[5])
#plt.plot(x_fit, y_fit, 'r',label='Fit Order: 5')
'''

# generate values for fit
x_fit = np.linspace(0,70,71)

# fit polynomial
polyfitting = np.polyfit(np.array(mapped_val), teff, 11)
poly11 = np.poly1d(polyfitting)

# PLOT
plt.clf()
fig, ax = plt.subplots()

# to adjust margin on left side of plot
fig.subplots_adjust(left=0.125)

# set limits
ax.set_ylim(-2,35000)
ax.set_xlim(-0.5,60)

# add some text for labels, title and axes ticks
ax.set_xlabel('Spectral Type')
ax.set_ylabel('$T_{eff}$')
### POSSIBLY: put titles inside plots, like rayner did
ax.set_title('Main Sequence: Astrophysical Quantities (Ed4)')


x_spec_list = ['B0', 'A0', 'F0', 'G0', 'K0', 'M0','S0']
ax.set_xticks([0,10,20,30,40,50,60])# + 1.5*width)
ax.set_xticklabels((x_spec_list))

# plot data
plt.plot(mapped_val, teff, 'ok', label='Lit. Values')
plt.plot(x_fit, poly11(x_fit), '--', label='Fit Order: 11')
#plt.legend()
print 'Teff of M0V: ', poly11(50)

# X-TICS
#ax = plt.axes()
#ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(2))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(1000))

plt.savefig('Interpolate_SptType_Teff.png')
#plt.show()
