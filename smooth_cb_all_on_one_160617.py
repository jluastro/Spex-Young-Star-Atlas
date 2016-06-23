# code3 mod of scripty_ intended to plot spectra of each star seperately
# Corey Mutnik
# 9/11/15-9/30/15 

import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.gridspec import GridSpec
from astropy.io import fits
from glob import glob
import matplotlib as mpl


#http://jdherman.github.io/colormap/


dir_path = u'/Users/cmutnik/work/astro/finished_with_fixed_names'
globpath = os.path.join(dir_path, '*.fits')
filelist = glob(globpath)
filelist.sort()


plt.clf()

nlines=len(filelist)

##
# Use the spectral colormap for examples
#cmap = plt.cm.Spectral

# Stacked_star_colors.png
cmap = mpl.colors.LinearSegmentedColormap.from_list('mycolors',['red','orange','yellow','white','blue'])

# Stacked_star_colors_nowhite.png
#cmap = mpl.colors.LinearSegmentedColormap.from_list('mycolors',['#ff0000', '#ff4500', '#ffa500', '#ffff00','#fffacd','#6495ed','#0000ff'])

# Use 0-1 values to generate the colors with the linspace method
line_colors = cmap(np.linspace(0,1,nlines))
colorbar_labels = np.linspace(1200, 50000, 21)
y = colorbar_labels


#fig = plt.figure(1)
# Generate fake ScalarMappable for colorbar
sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=y[0],vmax=y[-1]))
sm.set_array([])  # You have to set a dummy-array for this to work...
cbar = plt.colorbar(sm)#, cax=cbax1)
cbar.set_label('Teff')
cbar.set_ticks(y)
cbar.set_ticklabels(['{:4.0f}'.format(yi) for yi in y]) # Make 'em nicer-looking

fig = plt.figure(1)
plt.xlim(0.7,2.55)
plt.ylim(0,160)


for j in range(len(filelist)):
    spectra = fits.getdata(filelist[j])
    norm_wave = np.where(abs(spectra[0] - 2.20) == min(abs(spectra[0]-2.20)))[0][0]
    norm_den = (float)((spectra[0][norm_wave] * spectra[1][norm_wave])**(-1))
    norm_flux = []
    for i in range(0, len(spectra[0])):
        norm_flux.append(spectra[0][i] * spectra[1][i] * norm_den + j*3)

    file_name = os.path.basename(filelist[j])[:-5] #[:-5] prints file_name removing '.fits' from each
    #print file_name
    figure_name = file_name + '.png'

    y_pos = 1 + j
    #to_plot_regions = np.arange(0.0, 3.0 , 0.1)

    # what does label here do w/o a legend
    #plt.plot(spectra[0], norm_flux, c=line_colors[j],lw=3,label='{:3.1f}'.format(y[j]))
    plt.plot(spectra[0], norm_flux, c=line_colors[j],lw=0.5)


fig.subplots_adjust(left=0.125)
#plt.axis('tight')
plt.xlabel('Wavelength ($\mu$m)')
plt.ylabel('$\lambda f_\lambda / \lambda f_\lambda (2.2\mu m) + $ constant')

plt.savefig('Stacked_star_colors.pdf')
#plt.show()










