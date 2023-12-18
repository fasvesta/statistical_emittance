# SPS distribution generation from xpart examples 
import numpy as np
import xtrack as xt
import xpart as xp
import sys
sys.path.append("./../statisticalEmittance/")
from statisticalEmittance import *


bunch_intensity = 1e11
sigma_z = 22.5e-2
n_part = int(1e6)
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
statisticalOptics=r.measure_bunch_moments(particles)

assert np.isclose(statisticalOptics['nemitt_x'],nemitt_x,atol=5e-8)
assert np.isclose(statisticalOptics['nemitt_y'],nemitt_y,atol=5e-9)
