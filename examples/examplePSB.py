import numpy as np
from cpymad.madx import Madx
import xtrack as xt
import xpart as xp
import sys
sys.path.append("./../statisticalEmittance/")
from statisticalEmittance import *

mad = Madx()
mad.call('PSB/madx/psb_injection_example.madx')

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
epsn_x = []
epsn_y = []


for ii in range(10):
    print(ii)
    tracker.track(p_gaussian)
    bunch_moments=r.measure_bunch_moments(p_gaussian)
    epsn_x.append(bunch_moments['nemitt_x'])
    epsn_y.append(bunch_moments['nemitt_y'])

print('epsn_x = ',epsn_x)
print('epsn_y = ',epsn_y)
