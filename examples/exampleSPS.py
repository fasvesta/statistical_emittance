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
n_part = int(5e5)
nemitt_x = 2e-6
nemitt_y = 2.5e-6

filename = ('../../xtrack/test_data/sps_w_spacecharge/line_no_spacecharge_and_particle.json')
with open(filename, 'r') as fid:
    ddd = json.load(fid)
tracker = xt.Tracker(line=xt.Line.from_dict(ddd['line']))
part_ref = xp.Particles.from_dict(ddd['particle'])
tracker.line.particle_ref = part_ref

particles = xp.generate_matched_gaussian_bunch(
         num_particles=n_part, total_intensity_particles=bunch_intensity,
         nemitt_x=nemitt_x, nemitt_y=nemitt_y, sigma_z=sigma_z,
         particle_ref=part_ref,
         tracker=tracker)


r=statisticalEmittance(particles)
epsn_x = []
epsn_y = []
for ii in range(5):
    print(ii)
    tracker.track(particles)
    r.setInputDistribution(particles)
    epsn_x.append(r.getNormalizedEmittanceX())
    epsn_y.append(r.getNormalizedEmittanceY())

print('epsn_x = ',epsn_x)
print('epsn_y = ',epsn_y)