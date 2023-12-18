# SPS distribution generation from xpart examples 
import json
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
statisticalOptics=r.measure_bunch_moments(particles, 'coupling')
assert np.isclose(statisticalOptics['coupling'],0,atol=5e-5)
assert np.isclose(statisticalOptics['emitt_4d'],statisticalOptics['nemitt_x']/r.beta0/r.gamma0*statisticalOptics['nemitt_y']/r.beta0/r.gamma0,atol=1e-16)

klqsa=1.5e-2
skews=[i for i in line.element_dict if i[:4]=='lqsa' and len(i)<12]
for i in skews:
    line.element_dict[i].ksl[1]=klqsa

particles_coupling = xp.generate_matched_gaussian_bunch(
         num_particles=n_part, total_intensity_particles=bunch_intensity,
         nemitt_x=nemitt_x, nemitt_y=nemitt_y, sigma_z=sigma_z,
         line=line)

statisticalOptics_coupling=r.measure_bunch_moments(particles_coupling, 'coupling')
assert np.isclose(statisticalOptics_coupling['coupling'],2,atol=5e-1)
assert np.isclose(statisticalOptics['emitt_4d'],statisticalOptics_coupling['emitt_4d'],atol=1e-16)

klqsa=1e-2
skews=[i for i in line.element_dict if i[:4]=='lqsa' and len(i)<12]
for i in skews:
    line.element_dict[i].ksl[1]=klqsa

particles_coupling = xp.generate_matched_gaussian_bunch(
         num_particles=n_part, total_intensity_particles=bunch_intensity,
         nemitt_x=nemitt_x, nemitt_y=nemitt_y, sigma_z=sigma_z,
         line=line)

statisticalOptics_coupling=r.measure_bunch_moments(particles_coupling, 'coupling')
assert np.isclose(statisticalOptics_coupling['coupling'],0.2,atol=5e-2)
assert np.isclose(statisticalOptics['emitt_4d'],statisticalOptics_coupling['emitt_4d'],atol=1e-16)