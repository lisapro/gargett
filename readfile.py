'''
Created on 15. jul. 2016

@author: ELP
'''

from netCDF4 import Dataset
from math import sqrt,fabs
import numpy as np
#import numpy.ma as ma
from shutil import copyfile
copyfile('LL12zax.nc', 'LL12zax-kz.nc')

#nc_file_out = 'B3zax-kz.nc'
nc_file_out = 'LL12zax-kz.nc'

fl = Dataset(nc_file_out, mode='a')
temp = fl.variables['temp'][:]
salt = fl.variables['salt'][:]
zax = fl.variables['zax'][:]
time = fl.variables['time'][:]
latc = fl.variables['latc'][:]
lonc = fl.variables['lonc'][:]

kz = fl.createVariable("kz","f4",("time",'zax')) # time,zax,latc,lonc

def svan(s, t, po):
    r3500 = 1028.1063
    r4 = 4.8314E-4
    dr350 = 28.106331

    p=po/10.
    sr= sqrt(fabs(s))

    r1= ((((6.536332E-9*t-1.120083E-6)*t+1.001685E-4)*
        t-9.095290E-3)*t+6.793952E-2)*t-28.263737
    r2= (((5.3875E-9*t-8.2467E-7)*t+7.6438E-5)*t-4.0899E-3)*t+8.24493E-1
    r3= (-1.6546E-6*t+1.0227E-4)*t-5.72466E-3

    sig=(r4*s + r3*sr + r2)*s +r1

    v350p=1.0/r3500
    sva=-sig*v350p/(r3500+sig)
    sigma= sig + dr350

    if p != 0.0: 
        e = (9.1697E-10*t+2.0816E-8)*t-9.9348E-7
        bw = (5.2787E-8*t-6.12293E-6)*t+3.47718E-5
        b = bw + e*s

        d= 1.91075E-4
        c = (-1.6078E-6*t-1.0981E-5)*t+2.2838E-3
        aw = ((-5.77905E-7*t+1.16092E-4)*t+1.43713E-3)*t-0.1194975
        a = (d*sr + c)*s + aw

        b1 = (-5.3009E-4*t+1.6483E-2)*t+7.944E-2
        a1 = ((-6.1670E-5*t+1.09987E-2)*t-0.603459)*t+54.6746
        kw = (((-5.155288E-5*t+1.360477E-2)*t-2.327105)*t + 148.4206)*t-1930.06
        ko = (b1*sr + a1)*s + kw

        dk = (b*p+a)*p+ko
        k35 = (5.03217E-5*p+3.359406)*p+21582.27
        gam=p/k35
        pk=1.0-gam
        sva = sva * pk + (v350p+sva)*p*dk/(k35*(k35+dk))

        v350p= v350p*pk
        dr35p=gam/v350p
        dvan= sva/(v350p*(v350p+sva))
        sigma = dr350 + dr35p -dvan
        return sigma
    else:
        return sigma


density_temp = np.zeros(temp.shape[0])
kz_temp = np.zeros((temp.shape[0],temp.shape[1]))
dz = np.zeros(temp.shape[1])

dz[:]=1.

for i in range(time.shape[0]):    
    for j in range(zax.shape[0]):
        density_temp[j] = svan(salt[i,j], temp[i,j], zax[j])
    for j in range(zax.shape[0]):        
        if j < (zax.shape[0]-1):
            kz_temp[i,j] = 0.5E-6 /(sqrt(9.81/
                                    (1000.+(density_temp[j]+density_temp[j+1])/2.)
                          *max(0.0000001,(fabs(density_temp[j+1]-density_temp[j])/(dz[j])))))
        else : 
            kz_temp[i,j] = kz_temp[i,j-1]

kz[:,:] = kz_temp
fl.close()
print ('*** SUCCESS writing  ncfile ')
