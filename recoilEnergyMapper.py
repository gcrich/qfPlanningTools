# 
# quick approximate recoil calculations for neutron elastic scattering measurements
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1 import ImageGrid
from matplotlib import rc
from matplotlib import rcParams
import matplotlib.ticker as ticker
from matplotlib import gridspec
from matplotlib.font_manager import findfont, FontProperties
from matplotlib._cm import cubehelix
import matplotlib


rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})


rcParams.update({'figure.autolayout':True})
rcParams['font.family']='sans-serif'
rcParams['font.sans-serif']='Helvetica'
rcParams['legend.fontsize']=13
rcParams['axes.labelsize']=13
rcParams['axes.titlesize']=13
rcParams['xtick.labelsize']=11
rcParams['ytick.labelsize']=11


class masses(object):
    """Constants related to masses
    """
    electron = 511 # keV/c^2
    deuteron = 1.8756e+06  # keV /c^2
    neutron = 939565.0  # keV/c^2
    he3 = 2.809414e6  # keV/c^2
    iodine = 118210.76 * 1e3
    sodium = 21414.8342 * 1e3
    germanium = 67652. * 1e3
    xenon = 122299. * 1e3
    
    
    
def getMappableContours(contour):
    '''
    get a scalar mappable that allows for a colorbar that is continuous, even if only using discrete contour levels
    
    from https://stackoverflow.com/questions/44498631/continuous-colorbar-with-contour-levels
    '''
    norm = matplotlib.colors.Normalize(vmin=contour.vmin, 
                                   vmax = contour.vmax)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=contour.cmap)
    sm.set_array([])
    return sm
    
    
def recoilEnergy(scatterAngle, energyIncident, massTarget , massIncident = masses.neutron):
    """Calculate the recoil energy for an elastic scatter with specified parameters"""
    leadingCoeff = 2 * energyIncident * massIncident**2 / (massIncident + massTarget)**2
    
    bracketed = (massTarget/massIncident + np.sin(scatterAngle)**2 - 
                 np.cos(scatterAngle) * np.sqrt((massTarget/massIncident)**2 - 
                        np.sin(scatterAngle)**2))
    return leadingCoeff*bracketed
    



targetMass = masses.xenon
targetAbbreviation = 'xe'



En = [70, 200, 400, 580]
name = []
for e in En:
    name.append('$E_n = {}$ keV'.format(e))
recoilAngles = np.linspace(0., 180, 180)
smallAngles = recoilAngles[recoilAngles <= 30]

fig, [axFull, axZoom] = plt.subplots(2)
for e, title in zip(En,name):
    axFull.plot(recoilAngles, recoilEnergy(recoilAngles/360*2*np.pi, e, targetMass), label=title)
    axZoom.plot(smallAngles, recoilEnergy(smallAngles/360*2*np.pi, e, targetMass), label=title)
axFull.set_ylabel('Recoil energy (keVnr)')
axZoom.set_xlabel('Scattering angle (degrees)')
axZoom.set_xlim(0, 30)
axZoom.set_ylim()
plt.legend(loc='upper left')
plt.draw()



incidentEnergies = np.linspace(50, 700, 65)
gridX,gridY = np.meshgrid(recoilAngles, incidentEnergies)




recoilEnergies = [[recoilEnergy(recoilAngle/360*2*np.pi, incidentEnergy, targetMass) for recoilAngle in recoilAngles] for incidentEnergy in incidentEnergies]
fig, ax = plt.subplots(figsize=(8.5,5.25))
ims = ax.imshow(recoilEnergies, interpolation='None', cmap='cubehelix_r',
                origin='bottom')
ax.set_xlabel('Recoil angle (degrees)')
ax.set_ylabel('Incident neutron energy (keV)')
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='3%', pad=0.05, aspect='auto')
fig.colorbar(ims,label='Recoil energy (keVnr)', cax=cax)
plt.draw()




recoilEnergyGrid = recoilEnergy(gridX/360*2*np.pi, gridY, targetMass)

recoilContours = [0.1, 0.25, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 15.0, 20., 25., 30, 35]
fig, ax = plt.subplots(figsize=(8.5,5.25))
cplot = ax.contour(gridX, gridY, recoilEnergyGrid, levels=recoilContours, 
                   zorder=10, cmap='plasma_r')
plt.clabel(cplot)
ax.grid(True)
ax.xaxis.set_minor_locator(ticker.MultipleLocator(10))
plt.xlabel('Scattering angle (degrees)')
plt.ylabel('Incident neutron energy (keV)')
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size='3%', pad=0.05, aspect='auto')
#plt.colorbar(cplot,extend='both', label='Recoil energy (keVnr)', cax=cax)
smap = getMappableContours(cplot)
plt.colorbar(smap,
             label='Recoil energy (keVnr)', cax=cax)
plt.savefig('kinematics_' + targetAbbreviation + '.pdf')
plt.draw()


# make a zoomed contour plot
incidentEnergies_zoom = np.linspace(50, 200, 10)
recoilAngles_zoom = np.linspace(0., 120, 120)
gridX_zoom,gridY_zoom = np.meshgrid(recoilAngles_zoom, incidentEnergies_zoom)
recoilEnergyGrid_zoom = recoilEnergy(gridX_zoom/360*2*np.pi, gridY_zoom, targetMass)



recoilContours_zoom = [0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1, 1.25, 1.5, 2, 3, 4, 5]
cfig_zoom, ax_zoom = plt.subplots(figsize=(8.5,5.25))
cplot_zoom = ax_zoom.contour(gridX_zoom, gridY_zoom, recoilEnergyGrid_zoom, 
                              levels=recoilContours_zoom, zorder=10,
                              cmap='plasma_r')
plt.clabel(cplot_zoom)
ax_zoom.grid(True)
ax_zoom.xaxis.set_minor_locator(ticker.MultipleLocator(10))
plt.xlabel('Scattering angle (degrees)')
plt.ylabel('Incident neutron energy (keV)')
divider = make_axes_locatable(ax_zoom)
cax_zoom = divider.append_axes("right", size='3%', pad=0.05, aspect='auto')

smap = getMappableContours(cplot_zoom)
#plt.colorbar(cplot_zoom,extend='both', label='Recoil energy (keVnr)', cax=cax_zoom)
plt.colorbar(smap,
             label='Recoil energy (keVnr)', cax=cax_zoom)
plt.savefig('kinematics-zoom_' + targetAbbreviation +'.pdf')
plt.draw()





plt.show()