import numpy as np
from cpymad.madx import Madx
import time
import xtrack as xt
import xpart as xp
import sys
sys.path.append("./../statisticalEmittance/")
from statisticalEmittance import *

mad = Madx()
mad.call('../examples/PSB/madx/psb_injection_test.madx')

line= xt.Line.from_madx_sequence(mad.sequence['psb'])
line.particle_ref=xp.Particles(mass0=xp.PROTON_MASS_EV,
                               gamma0=mad.sequence.psb.beam.gamma)
tracker = xt.Tracker(line=line)

nemitt_x=3e-6
nemitt_y=1e-6
p_gaussian = xp.generate_matched_gaussian_bunch(num_particles=500000,
                            total_intensity_particles=5e11,
                            nemitt_x=nemitt_x, nemitt_y=nemitt_y, sigma_z=20.,
                            particle_ref=line.particle_ref,
                            tracker=tracker)



r=StatisticalEmittance()
t0=[]
epsn_x = []
epsn_y = []
for k in range(3):
    start = time.time()
    for ii in range(10):
        print(ii)
        tracker.track(p_gaussian)
        bunch_moments=r.measure_bunch_moments(p_gaussian)
        epsn_x.append(bunch_moments['nemitt_x'])
        epsn_y.append(bunch_moments['nemitt_y'])
    end = time.time()
    t0.append([end - start])
t1=[]
for k in range(3):
    start = time.time()
    for ii in range(10):
        print(ii)
        tracker.track(p_gaussian)
        bunch_moments=r.measure_bunch_moments(p_gaussian, coupling=True)
        epsn_x.append(bunch_moments['nemitt_x'])
        epsn_y.append(bunch_moments['nemitt_y'])
    end = time.time()
    t1.append([end - start])
t=[]
for k in range(3):
    start = time.time()
    for ii in range(10):
        print(ii)
        tracker.track(p_gaussian)
    end = time.time()
    t.append([end - start])

print('time tracking & estimating emittances (s)',t0)
print('time tracking & estimating emittances & coupling factor (s)',t1)
print('time tracking only (s)',t)
print('tracking & estimating emittances - only tracking (s)', np.array(t0)-np.array(t))