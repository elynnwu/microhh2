import matplotlib.pyplot as plt
from netCDF4 import Dataset
import  seaborn as sns
sns.set(style='ticks',font_scale=1.2)
xcross = 5
ycross = 5
f = Dataset('rflx.nc')
x = f['x'][:]
y = f['y'][:]
z = f['z'][:]
rflx = f['rflx']
plt.figure(figsize=(10,5))
plt.subplot(141)
plt.plot(rflx[0,:,xcross,ycross],z)
plt.ylim([0,2000])
plt.xlabel(r'Total radiative flux $[W/m^2]$')
g = Dataset('ql.nc')
plt.subplot(142)
ql = g['ql']
plt.plot(ql[0,:,xcross,ycross]*1000.,z)
plt.ylim([0,2000])
plt.xlabel(r'$q_l [g/kg]$')

f2 = Dataset('thl.nc')
thl = f2['thl']
plt.subplot(143)
plt.plot(thl[0,:,xcross,ycross],z)
plt.xlabel(r'$\theta_l [kg/kg]$')
plt.ylim([0,2000])

f3 = Dataset('qt.nc')
qt = f3['qt']
plt.subplot(144)
plt.plot(qt[0,:,xcross,ycross]*1000.,z)
plt.xlabel(r'$q_t [g/kg]$')
plt.ylim([0,2000])
plt.tight_layout()
# plt.savefig('dycoms_profile_at_x'+str(xcross)+'_y'+str(ycross)+'.png',dpi=200,bbox_inches='tight')
