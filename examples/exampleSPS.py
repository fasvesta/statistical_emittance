# SPS distribution generation from xpart examples 
import xtrack as xt
import xpart as xp
import sys
sys.path.append("./../statisticalEmittance/")
from statisticalEmittance import *


bunch_intensity = 1e11
sigma_z = 22.5e-2
n_part = int(1e5)
nemitt_x = 2e-6
nemitt_y = 2.5e-6

#filename = ('../../xtrack/test_data/sps_w_spacecharge/line_no_spacecharge_and_particle.json')
filename = ('line_no_spacecharge_and_particle.json')
line = xt.Line.from_json(filename)
line.build_tracker()


particles = xp.generate_matched_gaussian_bunch(
         num_particles=n_part, total_intensity_particles=bunch_intensity,
         nemitt_x=nemitt_x, nemitt_y=nemitt_y, sigma_z=sigma_z,
         line=line)


r=StatisticalEmittance()
epsn_x = []
epsn_y = []

for ii in range(10):
    print(ii)
    line.track(particles)
    bunch_moments=r.measure_bunch_moments(particles)
    epsn_x.append(bunch_moments['nemitt_x'])
    epsn_y.append(bunch_moments['nemitt_y'])

print('epsn_x = ',epsn_x)
print('epsn_y = ',epsn_y)