import numpy as np

from cpymad.madx import Madx

import time
import sys
sys.path.append("/home/fasvesta/statistical_emittance/statisticalEmittance")
from statisticalEmittance import *


mad = Madx()
mad.call('PSB/madx/psb_injection_for_pyOrbitNoErrors1.madx')


import xtrack as xt
import xpart as xp

line= xt.Line.from_madx_sequence(mad.sequence['psb'])
line.particle_ref=xp.Particles(mass0=xp.PROTON_MASS_EV,
                               gamma0=mad.sequence.psb.beam.gamma)
tracker = xt.Tracker(line=line)

tw = tracker.twiss()
mad.twiss()
beta0 = line.particle_ref.beta0[0]
print(f"Q'x (MAD)={mad.table.summ['dq1']*beta0} Q'x (Xsuite)={tw['dqx']}" )
print(f"Q'y (MAD)={mad.table.summ['dq2']*beta0} Q'y (Xsuite)={tw['dqy']}" )


p_gaussian = xp.generate_matched_gaussian_bunch(num_particles=500000,
                            total_intensity_particles=5e11,
                            nemitt_x=3e-6, nemitt_y=1e-6, sigma_z=20.,
                            particle_ref=line.particle_ref,
                            tracker=tracker)


g=p_gaussian.gamma0[0]
b=p_gaussian.beta0[0]

r=statisticalEmittance(x=p_gaussian.x,px=p_gaussian.px,y=p_gaussian.y,py=p_gaussian.py,z=p_gaussian.zeta,dp=p_gaussian.delta)
values=[]
values.append([r.getNormalizedEmittanceX(b,g),r.getNormalizedEmittanceY(b,g),r.getBetaX(),r.getBetaY()])

t=[]
for k in range(3):
    start = time.time()
    for ii in range(10):
        print(ii)
        tracker.track(p_gaussian)
        r=statisticalEmittance(x=p_gaussian.x,px=p_gaussian.px,y=p_gaussian.y,py=p_gaussian.py,z=p_gaussian.zeta,dp=p_gaussian.delta)
        values.append([r.getNormalizedEmittanceX(b,g),r.getNormalizedEmittanceY(b,g),r.getBetaX(),r.getBetaY()])
    end = time.time()
    t.append([end - start])
    np.save('values_3_2_'+str(k),values)

for k in range(3):
    start = time.time()
    for ii in range(10):
        print(ii)
        tracker.track(p_gaussian)
        #r=statisticalEmittance(x=p_gaussian.x,px=p_gaussian.px,y=p_gaussian.y,py=p_gaussian.py,z=p_gaussian.zeta,dp=p_gaussian.delta)
        #values.append([r.getNormalizedEmittanceX(b,g),r.getNormalizedEmittanceY(b,g),r.getBetaX(),r.getBetaY()])
    end = time.time()
    t.append([end - start])
    #np.save('values_3_2_'+str(k),values)
np.save('testTime',t)
