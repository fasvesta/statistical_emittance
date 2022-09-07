import numpy as np
import sys
sys.path.append("/home/fasvesta/statistical_emittance/statisticalEmittance")
from statisticalEmittance import *

# distribution generated from xpart examples 004_generate_gaussian.py
x=np.load('sps/datax.npy')
y=np.load('sps/datay.npy')
z=np.load('sps/dataz.npy')
dp=np.load('sps/datad.npy')
px=np.load('sps/datapx.npy')
py=np.load('sps/datapy.npy')
b=np.load('sps/datab.npy')[0]
g=np.load('sps/datag.npy')[0]

print('SPS input distribution values: NormalizedEmittanceX: 2e-6')
print('SPS input distribution values: NormalizedEmittanceY: 2.5e-6')

r=statisticalEmittance(x=x,px=px,y=y,py=py,z=z,dp=dp)
r.getFullOptics()
print('NormalizedEmittanceX: ', r.getNormalizedEmittanceX(b,g) )
print('NormalizedEmittanceY: ', r.getNormalizedEmittanceY(b,g) )
print('Optics',  r.getFullOptics())

# distribution generated from cpymadtest2.py

x=np.load('PSB/madx/datax.npy')
y=np.load('PSB/madx/datay.npy')
z=np.load('PSB/madx/dataz.npy')
dp=np.load('PSB/madx/datad.npy')
px=np.load('PSB/madx/datapx.npy')
py=np.load('PSB/madx/datapy.npy')
b=np.load('PSB/madx/datab.npy')[0]
g=np.load('PSB/madx/datag.npy')[0]

print('PSB input distribution values: NormalizedEmittanceX: 2e-6')
print('PSB input distribution values: NormalizedEmittanceY: 1e-6')
PSB_Optics= {'betx': 5.88, 'bety': 4.32, 'alfx': 2.41e-1, 'alfy': 3.58e-1, 'dispx': -2.81*b}
print('PSB madx optics: ', PSB_Optics)

r=statisticalEmittance(x=x,px=px,y=y,py=py,z=z,dp=dp)
r.getFullOptics()
print('NormalizedEmittanceX: ', r.getNormalizedEmittanceX(b,g) )
print('NormalizedEmittanceY: ', r.getNormalizedEmittanceY(b,g) )
print('Optics',  r.getFullOptics())

x=np.load('PSB/madx/datax3.npy')
y=np.load('PSB/madx/datay3.npy')
z=np.load('PSB/madx/dataz3.npy')
dp=np.load('PSB/madx/datad3.npy')
px=np.load('PSB/madx/datapx3.npy')
py=np.load('PSB/madx/datapy3.npy')


print('PSB input distribution values: NormalizedEmittanceX: 3e-6')
print('PSB input distribution values: NormalizedEmittanceY: 4e-6')
PSB_Optics= {'betx': 5.88, 'bety': 4.32, 'alfx': 2.41e-1, 'alfy': 3.58e-1, 'dispx': -2.81*b}
print('PSB madx optics: ', PSB_Optics)

r=statisticalEmittance(x=x,px=px,y=y,py=py,z=z,dp=dp)
r.getFullOptics()
print('NormalizedEmittanceX: ', r.getNormalizedEmittanceX(b,g) )
print('NormalizedEmittanceY: ', r.getNormalizedEmittanceY(b,g) )
print('Optics',  r.getFullOptics())
