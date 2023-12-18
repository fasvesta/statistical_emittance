from cpymad.madx import Madx

import xpart as xp
import xtrack as xt


seq_name = 'sps'
bunch_intensity = 1e11/3
sigma_z = 22.5e-2/3 # Short bunch to avoid probing bucket non-linearity
                    # to compare against frozen
nemitt_x=2.5e-6
nemitt_y=2.5e-6

mad = Madx()
mad.call('sps_thin.seq')
mad.use(seq_name)

line = xt.Line.from_madx_sequence(
                                            mad.sequence[seq_name],
                                            deferred_expressions=True,
                                            install_apertures=True)
# enable RF
V_RF = 3e6
line['acta.31637'].voltage = V_RF
line['acta.31637'].lag = 180.

# A test particle

line.particle_ref = xp.Particles(mass0=xp.PROTON_MASS_EV, gamma0=mad.sequence[seq_name].beam.gamma)

line.to_json('line_no_spacecharge_and_particle.json')
