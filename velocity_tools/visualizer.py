from astropy.constants import G
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord, FK5
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
import stream_lines as SL
from astropy.io import fits
from matplotlib.widgets import Slider, Button
from scipy import stats
from regions import read_ds9
import pyregion


class Visualizer():

    image = None
    imageheader = None
    star_c = None
    star_ref = None
    zz = None

    def __init__(self, v_lsr=7.5 *u.km/u.s, Mstar=1 * u.Msun, omega0=1e-13 / u.s, theta0=45 * u.deg,
                 phi0=0 * u.deg, v_r0=0 * u.km/u.s, r0=1000 * u.au,
                 inc=0 * u.deg, PA_ang =0 * u.deg, rmin_stream = 100*u.au):
        # initial params of the model
        self.v_lsr = v_lsr
        self.Mstar = Mstar
        self.omega0 = omega0
        self.theta0 = theta0
        self.phi0 = phi0
        self.v_r0 = v_r0
        self.r0 = r0
        self.inc = inc
        self.PA_ang = PA_ang
        self.rmin_stream = rmin_stream

        self.fig = plt.figure(figsize=(10,7))
        plt.subplots_adjust(left=0.1, bottom=0.45)

    def set_referenceframe(self, centerra, centerdec, distance):
        self.star_c = SkyCoord(centerra, centerdec, frame='fk5')
        self.star_ref = self.star_c.skyoffset_frame()
        self.distance = distance

    def set_image(self, infile, vmin, vmax):
        if isinstance(infile, str):
            self.imagedata = fits.getdata(infile)
            self.imageheader = fits.getheader(infile)
        else:
            self.imagedata = infile.data
            self.imageheader = infile.header
        self.vmin = vmin
        self.vmax = vmax

    def get_vc_r(self, velfield_file, region_file):
        """
        Returns the centroid velocity and projected separation in the sky for the
        centroid velocity from any star

        Given a region and a velocity field for the vicinity of a star,
        obtains the projected radius and the central velocity of each pixel in the
        region. The velocity field must be masked to contain only the relevant
        pixels.

        Args:
            velfield_file (string): path to the .fits file containing the velocity
            field
            region_file (string): path to the .reg (ds9) region file where the
            streamer is contained

        Returns:
            type: description

        """
        import velocity_tools.coordinate_offsets as c_offset
        # load region file and WCS structures
        regions = read_ds9(region_file)
        wcs_Vc = WCS(velfield_file)
        #
        hd_Vc = fits.getheader(velfield_file)
        results = c_offset.generate_offsets(
            hd_Vc, self.star_c.ra, self.star_c.dec, pa_angle=0*u.deg, inclination=0*u.deg)
        rad_au = (results.r * self.distance*u.pc).to(u.au, equivalencies=u.dimensionless_angles())

        mask_Vc = (regions[0].to_pixel(wcs_Vc)).to_mask()
        Vc_cutout = mask_Vc.cutout(fits.getdata(velfield_file))
        rad_cutout = mask_Vc.cutout(rad_au)
        #
        gd = (mask_Vc.data == 1)
        v_los = Vc_cutout[gd]*u.km/u.s
        r_proj = rad_cutout[gd]
        return r_proj, v_los


    def set_velocity_kde(self, velimagename, regionsamplepath, rmin=0, rmax=1000, ymin=None, ymax=None):
        self.regionname = regionsamplepath
        r_proj, v_los = self.get_vc_r(velimagename, self.regionname)
        # x is projected distance in pc
        # y is velocity lsr in kms
        self.velmin = self.v_lsr.value-2 if ymin is None else ymin
        self.velmax = self.v_lsr.value+2 if ymax is None else ymax
        self.rmin = rmin
        self.rmax = rmax

        self.xx, self.yy = np.mgrid[self.rmin:self.rmax:100j, self.velmin:self.velmax:100j]
        positions = np.vstack([self.xx.ravel(), self.yy.ravel()])
        gd_vlos = np.isfinite(r_proj*v_los)
        values = np.vstack([r_proj[gd_vlos].value, v_los[gd_vlos].value])
        # we calculate the kernel distribution
        kernel = stats.gaussian_kde(values)
        self.zz = np.reshape(kernel(positions).T, self.xx.shape)
        self.zz /= self.zz.max()  # normalization of probability
        # xx, yy and zz are our kde.



    def calculate_streamline(self):
        (x1, y1, z1), (vx1, vy1, vz1) = SL.xyz_stream(mass=self.Mstar, r0=self.r0,
                                                      theta0=self.theta0, phi0=self.phi0,
                                                      omega=self.omega0, v_r0=self.v_r0,
                                                      inc=self.inc, pa=self.PA_ang,
                                                      rmin=self.rmin_stream)
        # we save the spatial and velocity info
        self.x1 = x1
        self.y1 = y1
        self.z1 = z1
        self.vx1 = vx1
        self.vy1 = vy1
        self.vz1 = vz1
        self.r_c = SL.r_cent(self.Mstar,self.omega0,self.r0)

    def calculate_image_plane_sl(self):
        '''
        Distance must be in pc as au/pc = arcsec
        '''
        dra_stream = -self.x1.value / self.distance
        ddec_stream = self.z1.value / self.distance
        self.d_sky_au = np.sqrt(self.x1**2 + self.z1**2)
        self.fil = SkyCoord(dra_stream*u.arcsec, ddec_stream*u.arcsec,
                       frame=self.star_ref).transform_to(FK5)
        self.velocity = self.v_lsr + self.vy1

    # TODO: do all the gets and sets

    def plot_sl(self, plotregion=True):
        wcs = WCS(self.imageheader)

        #ax is the image plane
        self.ax = self.fig.add_subplot(121, projection=wcs)
        self.ax.imshow(self.imagedata, vmin=self.vmin, vmax=self.vmax, origin='lower', cmap='Greys')
        self.ax.set_xlabel('Right Ascension (J2000)')
        self.ax.set_ylabel('Declination (J2000)')
        if plotregion:
            regstreamer = pyregion.open(self.regionname)
            r2 = regstreamer.as_imagecoord(self.imageheader)
            patch_list, artist_list = r2.get_mpl_patches_texts()
            for p in patch_list:
                self.ax.add_patch(p)
            for a in artist_list:
                self.ax.add_artist(a)

        # ax2 is the velocity kde plane
        self.ax2 = self.fig.add_subplot(122)
        self.ax2.set_xlabel('Projected distance (au)')
        self.ax2.set_ylabel(r"V$_{lsr}$ (km s$^{-1}$)")
        if self.zz is not None:
            self.ax2.contourf(self.xx, self.yy, self.zz, cmap='Greys', levels=np.arange(0.1, 1.2, 0.1), vmin=0., vmax=1.1)
            self.ax2.set_ylim([self.velmin, self.velmax])
            self.ax2.set_xlim([self.rmin, self.rmax])

        # here we plot the streamline model
        self.line_image, = self.ax.plot(self.fil.ra, self.fil.dec, transform=self.ax.get_transform('fk5'),
                      ls='-', lw=2)
        self.line_vel, = self.ax2.plot(self.d_sky_au, self.velocity)
        self.annotation = self.ax2.annotate(r'$r_c = {}$'.format(np.round(self.r_c,0)), (0.6, 0.1), xycoords='axes fraction', size=12)

        delta_theta0 = 5
        delta_phi0 = 5
        delta_r0 = 10
        delta_omega0 = 1.e-14
        delta_v_r0 = 0.1

        axcolor = 'paleturquoise'

        axtheta0 = plt.axes([0.2, 0.1, 0.6, 0.03], facecolor=axcolor)
        self.stheta0 = Slider(axtheta0, r'$\theta_0$', 0, 180., valinit=self.theta0.value, valstep=delta_theta0)

        axphi0 = plt.axes([0.2, 0.15, 0.6, 0.03], facecolor=axcolor)
        self.sphi0 = Slider(axphi0, r'$\phi_0$', 0, 360, valinit=self.phi0.value, valstep=delta_phi0)

        axr0 = plt.axes([0.2, 0.2, 0.6, 0.03], facecolor=axcolor)
        self.sr0 = Slider(axr0, r'$r_0$', 1500, 10000, valinit=self.r0.value, valstep=delta_r0)

        axomega0 = plt.axes([0.2, 0.25, 0.6, 0.03], facecolor=axcolor)
        self.somega0 = Slider(axomega0, r'$\Omega_0$', 1.e-14, 1e-12, valinit=self.omega0.value, valstep=delta_omega0)

        axv0 = plt.axes([0.2, 0.3, 0.6, 0.03], facecolor=axcolor)
        self.sv0 = Slider(axv0, r'$v_{r,0}$', 0, 5, valinit=self.v_r0.value, valstep=delta_v_r0)

        updateax = plt.axes([0.8, 0.025, 0.1, 0.04])
        button = Button(updateax, 'Update', color=axcolor, hovercolor='0.975')
        button.on_clicked(self.update)

        plt.show()

    def update(self, val):
        self.theta0 = self.stheta0.val * u.deg
        self.phi0 = self.sphi0.val * u.deg
        self.omega0 = self.somega0.val / u.s
        self.r0 = self.sr0.val * u.au
        self.v_r0 = self.sv0.val * u.km / u.s
        self.calculate_streamline()
        self.calculate_image_plane_sl()

        self.line_image.set_xdata(self.fil.ra)
        self.line_image.set_ydata(self.fil.dec)
        self.line_vel.set_xdata(self.d_sky_au)
        self.line_vel.set_ydata(self.velocity)
        self.r_c = SL.r_cent(self.Mstar, self.omega0, self.r0)
        self.annotation.set_text(r'$r_c = {}$'.format(np.round(self.r_c,0)))
        self.fig.canvas.draw_idle()
