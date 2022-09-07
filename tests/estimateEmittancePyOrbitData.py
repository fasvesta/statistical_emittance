import numpy as np
import scipy.io as sio 
import sys
sys.path.append("/home/fasvesta/statistical_emittance/statisticalEmittance")
from statisticalEmittance import *

# distributions from old simulation data /eos/project/l/liu/PSB/spaceChargeSimulationData/2018/tuneScan2QyOn/mdTuneScans4p524
data=np.loadtxt('PSB/madx/mainbunch_start1.dat',comments='%')
b=np.load('PSB/madx/datab.npy')[0]
g=np.load('PSB/madx/datag.npy')[0]

print('pyorbit : NormalizedEmittanceX: 3.01385982021e-06')
print('pyorbit : NormalizedEmittanceY: 2.29405965534e-06')
PSB_Optics= {'betx': 6.11, 'bety': 4.22, 'alfx': 2.31e-1, 'alfy': 3.68e-1, 'dispx': -1.49}
print('PSB optics: ', PSB_Optics)

r=statisticalEmittance(x=data[:,0],px=data[:,1],y=data[:,2],py=data[:,3],z=data[:,4],dp=data[:,5]) 

print('NormalizedEmittanceX: ', r.getNormalizedEmittanceX(b,g) )
print('NormalizedEmittanceY: ', r.getNormalizedEmittanceY(b,g) )
print('Optics ',  r.getFullOptics())
print('corrected dispersion: ', r.dispersionX*0.5708*b)

data=sio.loadmat('PSB/madx/mainbunch_002999.mat') 

print('pyorbit : NormalizedEmittanceX: 3.01536537217e-06')
print('pyorbit : NormalizedEmittanceY: 2.86442324125e-06')

r=statisticalEmittance(x=data['particles'][0][0]['x'].flatten(),px=data['particles'][0][0]['xp'].flatten(),y=data['particles'][0][0]['y'].flatten(),py=data['particles'][0][0]['yp'].flatten(),z=data['particles'][0][0]['z'].flatten(),dp=data['particles'][0][0]['dE'].flatten())

print('NormalizedEmittanceX: ', r.getNormalizedEmittanceX(b,g) )
print('NormalizedEmittanceY: ', r.getNormalizedEmittanceY(b,g) )
print('Optics',  r.getFullOptics())

data=sio.loadmat('PSB/madx/mainbunch_000999.mat')

r=statisticalEmittance(x=data['particles'][0][0]['x'].flatten(),px=data['particles'][0][0]['xp'].flatten(),y=data['particles'][0][0]['y'].flatten(),py=data['particles'][0][0]['yp'].flatten(),z=data['particles'][0][0]['z'].flatten(),dp=data['particles'][0][0]['dE'].flatten())

print('pyorbit : NormalizedEmittanceX: 3.01436715678e-06')
print('pyorbit : NormalizedEmittanceY: 2.51335244655e-06')

print('NormalizedEmittanceX: ', r.getNormalizedEmittanceX(b,g) )
print('NormalizedEmittanceY: ', r.getNormalizedEmittanceY(b,g) )
print('Optics',  r.getFullOptics())

data=sio.loadmat('PSB/q4q3/mainbunch_009999.mat')

r=statisticalEmittance(x=data['particles'][0][0]['x'].flatten(),px=data['particles'][0][0]['xp'].flatten(),y=data['particles'][0][0]['y'].flatten(),py=data['particles'][0][0]['yp'].flatten(),z=data['particles'][0][0]['z'].flatten(),dp=data['particles'][0][0]['dE'].flatten())

data1=sio.loadmat('PSB/q4q3/output.mat')
print('pyorbit : NormalizedEmittanceX: ', data1['epsn_x'][0][-1])
print('pyorbit : NormalizedEmittanceY: ', data1['epsn_y'][0][-1])

print('NormalizedEmittanceX: ', r.getNormalizedEmittanceX(b,g) )
print('NormalizedEmittanceY: ', r.getNormalizedEmittanceY(b,g) )
print('Optics',  r.getFullOptics())
print('Coupling',  r.getCouplingFactor())

data=sio.loadmat('PSB/q4q4/mainbunch_009999.mat')

r=statisticalEmittance(x=data['particles'][0][0]['x'].flatten(),px=data['particles'][0][0]['xp'].flatten(),y=data['particles'][0][0]['y'].flatten(),py=data['particles'][0][0]['yp'].flatten(),z=data['particles'][0][0]['z'].flatten(),dp=data['particles'][0][0]['dE'].flatten())

data1=sio.loadmat('PSB/q4q4/output.mat')
print('pyorbit : NormalizedEmittanceX: ', data1['epsn_x'][0][-1])
print('pyorbit : NormalizedEmittanceY: ', data1['epsn_y'][0][-1])

print('NormalizedEmittanceX: ', r.getNormalizedEmittanceX(b,g) )
print('NormalizedEmittanceY: ', r.getNormalizedEmittanceY(b,g) )
print('Optics',  r.getFullOptics())
print('Coupling',  r.getCouplingFactor())

data=sio.loadmat('PS/mainbunch_499999.mat')
b=0.91596101476611
g=2.4921045962752
p0=2.1417660241115 

r=statisticalEmittance(x=data['particles'][0][0]['x'].flatten(),px=data['particles'][0][0]['xp'].flatten(),y=data['particles'][0][0]['y'].flatten(),py=data['particles'][0][0]['yp'].flatten(),z=data['particles'][0][0]['z'].flatten(),dp=data['particles'][0][0]['dE'].flatten())

data1=sio.loadmat('PS/output_-1.mat')
print('pyorbit : NormalizedEmittanceX: ', data1['epsn_x'][0][-1])
print('pyorbit : NormalizedEmittanceY: ', data1['epsn_y'][0][-1])
PS_Optics={'betx': 21.58, 'bety': 11.78, 'alfx': 4.53e-2, 'alfy': 0.1, 'dispx': 3.41}
print('PS optics: ', PS_Optics)

print('NormalizedEmittanceX: ', r.getNormalizedEmittanceX(b,g) )
print('NormalizedEmittanceY: ', r.getNormalizedEmittanceY(b,g) )
print('Optics',  r.getFullOptics())
print('corrected dispersion: ', r.dispersionX*p0*b)
