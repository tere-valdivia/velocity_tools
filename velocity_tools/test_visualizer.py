from visualizer import *

fitsfile = '../examples/JEP_testdata/Per-emb-2-HC3N_10-9_TdV'
velimage = '../examples/JEP_testdata/Per-emb-2-HC3N_10-9_fit_Vc'
regionsamplepath = '../examples/JEP_testdata/Streamer_North_v2.reg'
ra_Per2 = 15 * (3 + (32 + 17.92/60.) / 60.) * u.deg
dec_Per2 = (30 + (49 + 48.03 / 60.) / 60.) * u.deg
distance = 300.  # pc

Mstar = 3.2*u.Msun
inc = -43*u.deg
PA_ang = 130*u.deg
theta0 = 130.*u.deg
r0 = 0.9e4*u.au
phi0 = 365.*u.deg
v_r0 = 0*u.km/u.s
omega0 = 4e-13/u.s
v_lsr = 7.05*u.km/u.s

viz = Visualizer(v_lsr=v_lsr, Mstar=Mstar, omega0=omega0, theta0=theta0, phi0=phi0, v_r0=v_r0, r0=r0, inc=inc, PA_ang = PA_ang, rmin_stream=5.5e3*u.au)
viz.set_image(fitsfile+'.fits', vmin=0, vmax=0.17)
viz.set_referenceframe(ra_Per2, dec_Per2, distance)
viz.set_velocity_kde(velimage+'.fits', regionsamplepath, rmax=9500)
viz.calculate_streamline()
viz.calculate_image_plane_sl()
viz.plot_sl()
'''

fitsfile = '../examples/JEP_testdata/Per-emb-2-HC3N_10-9_TdV'
velimage = '../examples/JEP_testdata/Per-emb-2-HC3N_10-9_fit_Vc'
regionsamplepath = '../examples/JEP_testdata/Streamer_North_v2.reg'
ra_Per2 = 15 * (3 + (32 + 17.92/60.) / 60.) * u.deg
dec_Per2 = (30 + (49 + 48.03 / 60.) / 60.) * u.deg
distance = 300.  # pc

Mstar = 3.2*u.Msun
inc = -43*u.deg
PA_ang = 130*u.deg
theta0 = 130.*u.deg
r0 = 9.e4*u.au
phi0 = 365.*u.deg
v_r0 = 0*u.km/u.s
omega0 = 4e-13/u.s
v_lsr = 7.05*u.km/u.s

viz = Visualizer(v_lsr=v_lsr, Mstar=Mstar, omega0=omega0, theta0=theta0, phi0=phi0, v_r0=v_r0, r0=r0, inc=inc, PA_ang = PA_ang)

'''
