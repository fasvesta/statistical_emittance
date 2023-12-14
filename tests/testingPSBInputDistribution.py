import numpy as np
from cpymad.madx import Madx
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

mad.twiss()
madOptics={'betx': mad.table.twiss['betx'][0],
 'bety': mad.table.twiss['bety'][0],
 'alfx': mad.table.twiss['alfx'][0],
 'alfy': mad.table.twiss['alfy'][0],
 'dx': mad.table.twiss['dx'][0],
 'dy': mad.table.twiss['dy'][0],
 'dpx': mad.table.twiss['dpx'][0],
 'dpy':mad.table.twiss['dpy'][0]}

nemitt_x=3e-6
nemitt_y=1e-6
p_gaussian = xp.generate_matched_gaussian_bunch(num_particles=500000,
                            total_intensity_particles=5e11,
                            nemitt_x=nemitt_x, nemitt_y=nemitt_y, sigma_z=20.,
                            particle_ref=line.particle_ref,
                            tracker=tracker)



r=StatisticalEmittance()
statisticalOptics=r.measure_bunch_moments(p_gaussian)

assert np.isclose(statisticalOptics['betx'],madOptics['betx'],atol=5e-2, rtol=0)
assert np.isclose(statisticalOptics['betx'],trackerOptics['betx'],atol=5e-2, rtol=0)

assert np.isclose(statisticalOptics['alfx'],madOptics['alfx'],atol=5e-3, rtol=0)
assert np.isclose(statisticalOptics['alfx'],trackerOptics['alfx'],atol=5e-3)

assert np.isclose(statisticalOptics['gamx'],trackerOptics['gamx'],atol=5e-4)

assert np.isclose(statisticalOptics['dx'],madOptics['dx']*p_gaussian.beta0[0],atol=5e-2)
assert np.isclose(statisticalOptics['dx'],trackerOptics['dx'],atol=5e-2)

assert np.isclose(statisticalOptics['dpx'],madOptics['dpx']*p_gaussian.beta0[0],atol=5e-3)
assert np.isclose(statisticalOptics['dpx'],trackerOptics['dpx'],atol=5e-3)

assert np.isclose(statisticalOptics['bety'],madOptics['bety'],atol=5e-2)
assert np.isclose(statisticalOptics['bety'],trackerOptics['bety'],atol=5e-2)

assert np.isclose(statisticalOptics['alfy'],madOptics['alfy'],atol=5e-3)
assert np.isclose(statisticalOptics['alfy'],trackerOptics['alfy'],atol=5e-3)

assert np.isclose(statisticalOptics['gamy'],trackerOptics['gamy'],atol=5e-3)

assert np.isclose(statisticalOptics['dy'],madOptics['dy']*p_gaussian.beta0[0],atol=5e-3)
assert np.isclose(statisticalOptics['dy'],trackerOptics['dy'],atol=5e-3)

assert np.isclose(statisticalOptics['dpy'],madOptics['dpy']*p_gaussian.beta0[0],atol=5e-3)
assert np.isclose(statisticalOptics['dpy'],trackerOptics['dpy'],atol=5e-3)

assert np.isclose(statisticalOptics['nemitt_x'],nemitt_x,atol=5e-8)
assert np.isclose(statisticalOptics['nemitt_y'],nemitt_y,atol=5e-9)


for ii in range(3):
    print(ii)
    tracker.track(p_gaussian)
    statisticalOptics=r.measure_bunch_moments(p_gaussian)
    assert np.isclose(statisticalOptics['betx'],madOptics['betx'],atol=5e-2)
    assert np.isclose(statisticalOptics['betx'],trackerOptics['betx'],atol=5e-2)

    assert np.isclose(statisticalOptics['alfx'],madOptics['alfx'],atol=5e-3)
    assert np.isclose(statisticalOptics['alfx'],trackerOptics['alfx'],atol=5e-3)

    assert np.isclose(statisticalOptics['gamx'],trackerOptics['gamx'],atol=5e-4)

    assert np.isclose(statisticalOptics['dx'],madOptics['dx']*p_gaussian.beta0[0],atol=5e-2)
    assert np.isclose(statisticalOptics['dx'],trackerOptics['dx'],atol=5e-2)

    assert np.isclose(statisticalOptics['dpx'],madOptics['dpx']*p_gaussian.beta0[0],atol=5e-3)
    assert np.isclose(statisticalOptics['dpx'],trackerOptics['dpx'],atol=5e-3)

    assert np.isclose(statisticalOptics['bety'],madOptics['bety'],atol=5e-2)
    assert np.isclose(statisticalOptics['bety'],trackerOptics['bety'],atol=5e-2)

    assert np.isclose(statisticalOptics['alfy'],madOptics['alfy'],atol=5e-3)
    assert np.isclose(statisticalOptics['alfy'],trackerOptics['alfy'],atol=5e-3)

    assert np.isclose(statisticalOptics['gamy'],trackerOptics['gamy'],atol=5e-3)

    assert np.isclose(statisticalOptics['dy'],madOptics['dy']*p_gaussian.beta0[0],atol=5e-3)
    assert np.isclose(statisticalOptics['dy'],trackerOptics['dy'],atol=5e-3)

    assert np.isclose(statisticalOptics['dpy'],madOptics['dpy']*p_gaussian.beta0[0],atol=5e-3)
    assert np.isclose(statisticalOptics['dpy'],trackerOptics['dpy'],atol=5e-3)

    assert np.isclose(statisticalOptics['nemitt_x'],nemitt_x,atol=5e-8)
    assert np.isclose(statisticalOptics['nemitt_y'],nemitt_y,atol=5e-9)