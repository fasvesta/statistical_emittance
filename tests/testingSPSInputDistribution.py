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

tw = tracker.twiss()
trackerOptics={'betx': tw['betx'][0],
 'bety': tw['bety'][0],
 'alfx': tw['alfx'][0],
 'alfy': tw['alfy'][0],
 'gamx': tw['gamx'][0],
 'gamy': tw['gamy'][0],
 'dx': tw['dx'][0],
 'dy': tw['dy'][0],
 'dpx': tw['dpx'][0],
 'dpy':tw['dpy'][0]}

particles = xp.generate_matched_gaussian_bunch(
         num_particles=n_part, total_intensity_particles=bunch_intensity,
         nemitt_x=nemitt_x, nemitt_y=nemitt_y, sigma_z=sigma_z,
         particle_ref=part_ref,
         tracker=tracker)


r=StatisticalEmittance(particles)
statisticalOptics=r.get_full_optics()

assert np.isclose(statisticalOptics['betx'],trackerOptics['betx'],atol=5e-1)

assert np.isclose(statisticalOptics['alfx'],trackerOptics['alfx'],atol=5e-3)

assert np.isclose(statisticalOptics['gamx'],trackerOptics['gamx'],atol=5e-4)

assert np.isclose(statisticalOptics['dx'],trackerOptics['dx'],atol=5e-2)

assert np.isclose(statisticalOptics['dpx'],trackerOptics['dpx'],atol=5e-3)

assert np.isclose(statisticalOptics['bety'],trackerOptics['bety'],atol=5e-1)

assert np.isclose(statisticalOptics['alfy'],trackerOptics['alfy'],atol=5e-3)

assert np.isclose(statisticalOptics['gamy'],trackerOptics['gamy'],atol=5e-3)

assert np.isclose(statisticalOptics['dy'],trackerOptics['dy'],atol=5e-3)

assert np.isclose(statisticalOptics['dpy'],trackerOptics['dpy'],atol=5e-3)

assert np.isclose(r.get_nemitt_x(),nemitt_x,atol=5e-8)
assert np.isclose(r.get_nemitt_y(),nemitt_y,atol=5e-9)


for ii in range(3):
    print(ii)
    tracker.track(particles)
    r.set_particles(particles)
    statisticalOptics=r.get_full_optics()
    assert np.isclose(statisticalOptics['betx'],trackerOptics['betx'],atol=5e-1)

    assert np.isclose(statisticalOptics['alfx'],trackerOptics['alfx'],atol=5e-2)

    assert np.isclose(statisticalOptics['gamx'],trackerOptics['gamx'],atol=5e-4)

    assert np.isclose(statisticalOptics['dx'],trackerOptics['dx'],atol=5e-2)

    assert np.isclose(statisticalOptics['dpx'],trackerOptics['dpx'],atol=5e-3)

    assert np.isclose(statisticalOptics['bety'],trackerOptics['bety'],atol=5e-1)

    assert np.isclose(statisticalOptics['alfy'],trackerOptics['alfy'],atol=5e-2)

    assert np.isclose(statisticalOptics['gamy'],trackerOptics['gamy'],atol=5e-3)

    assert np.isclose(statisticalOptics['dy'],trackerOptics['dy'],atol=5e-3)

    assert np.isclose(statisticalOptics['dpy'],trackerOptics['dpy'],atol=5e-3)

    assert np.isclose(r.get_nemitt_x(),nemitt_x,atol=5e-8)
    assert np.isclose(r.get_nemitt_y(),nemitt_y,atol=5e-9)