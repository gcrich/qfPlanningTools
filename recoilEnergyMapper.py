# 
# quick approximate recoil calculations for neutron elastic scattering measurements
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt


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
    
    
def recoilEnergy(scatterAngle, energyIncident, massTarget , massIncident = masses.neutron):
    """Calculate the recoil energy for an elastic scatter with specified parameters"""
    leadingCoeff = 2 * energyIncident * massIncident**2 / (massIncident + massTarget)**2
    
    bracketed = (massTarget/massIncident + np.sin(scatterAngle)**2 - 
                 np.cos(scatterAngle) * np.sqrt((massTarget/massIncident)**2 - 
                        np.sin(scatterAngle)**2))
    return leadingCoeff*bracketed
    
En = [70, 200, 400, 580]
name = []
for e in En:
    name.append('$E_n = {}$ keV'.format(e))
recoilAngles = np.linspace(0., 180, 180)
smallAngles = recoilAngles[recoilAngles <= 30]

fig, [axFull, axZoom] = plt.subplots(2)
for e, title in zip(En,name):
    axFull.plot(recoilAngles, recoilEnergy(recoilAngles/360*2*np.pi, e, masses.germanium), label=title)
    axZoom.plot(smallAngles, recoilEnergy(smallAngles/360*2*np.pi, e, masses.germanium), label=title)
axFull.set_ylabel('Recoil energy (keVnr)')
axZoom.set_xlabel('Scattering angle (degrees)')
axZoom.set_xlim(0, 30)
axZoom.set_ylim()
plt.legend(loc='upper left')
plt.draw()



incidentEnergies = np.linspace(50, 700, 65)
gridX,gridY = np.meshgrid(recoilAngles, incidentEnergies)
recoilEnergies = [[recoilEnergy(recoilAngle/360*2*np.pi, incidentEnergy, masses.germanium) for recoilAngle in recoilAngles] for incidentEnergy in incidentEnergies]
fig, ax = plt.subplots()
ims = ax.imshow(recoilEnergies, interpolation='None', cmap='viridis')
fig.colorbar(ims,ax=ax)
plt.draw()

recoilEnergyGrid = recoilEnergy(gridX/360*2*np.pi, gridY, masses.germanium)

recoilContours = [0.1, 0.25, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 15.0, 20., 25., 30, 35]
fig, ax = plt.subplots()
cplot = ax.contour(gridX, gridY, recoilEnergyGrid, levels=recoilContours)
plt.clabel(cplot)
ax.grid(True)
plt.xlabel('Scattering angle (degrees)')
plt.ylabel('Incident neutron energy (keV)')
plt.colorbar(cplot,extend='both', label='Recoil energy (keVnr)')
#plt.savefig('ge-kinematics.png',dpi=400)
plt.draw()


# make a zoomed contour plot
incidentEnergies_zoom = np.linspace(50, 200, 10)
recoilAngles_zoom = np.linspace(0., 120, 120)
gridX_zoom,gridY_zoom = np.meshgrid(recoilAngles_zoom, incidentEnergies_zoom)
recoilEnergyGrid_zoom = recoilEnergy(gridX_zoom/360*2*np.pi, gridY_zoom, masses.germanium)
recoilContours_zoom = [0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1, 1.25, 1.5, 2, 3, 4]
cfig_zoom, cax_zoom = plt.subplots()
cplot_zoom = cax_zoom.contour(gridX_zoom, gridY_zoom, recoilEnergyGrid_zoom, levels=recoilContours_zoom)
plt.clabel(cplot_zoom)
cax_zoom.grid(True)
plt.xlabel('Scattering angle (degrees)')
plt.ylabel('Incident neutron energy (keV)')
plt.colorbar(cplot_zoom,extend='both', label='Recoil energy (keVnr)')
#plt.savefig('ge-kinematics-zoom.png',dpi=400)
plt.draw()





plt.show()