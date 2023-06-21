import os
import glob
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)
import numpy as np
import math
import matplotlib
matplotlib.use('Agg')
import sys
from jet_image import JetImage, TwinJetImage
from vlbi_utils import find_image_std, find_bbox, pol_mask, correct_ppol_bias
sys.path.insert(0, '/home/ilya/github/ve/vlbi_errors')
from uv_data import UVData
from spydiff import clean_difmap, find_nw_beam
from from_fits import create_clean_image_from_fits_file
from image import plot as iplot
sys.path.insert(0, '/home/ilya/github/stackemall')
from stack_utils import stat_of_masked
import matplotlib.pyplot as plt


def convert_mojave_epoch_to_float(epoch):
    year = epoch.split('_')[0]
    month = epoch.split('_')[1]
    day = epoch.split('_')[2]
    result = float(year) + float(month)/12.0 + float(day)/30.0
    return result


def get_epochs_from_image_list(files):
    times = list()
    epochs = list()
    for file in files:
        file = os.path.split(file)[-1]
        epoch = file.split(".")[2]
        epochs.append(epoch)
        time = convert_mojave_epoch_to_float(epoch)
        times.append(time)
    times = np.array(times)
    times = np.round(times - np.min(times), 2)
    return times, epochs


def plot_function(contours=None, colors=None, vectors=None, vectors_values=None,
                  cmap='gist_rainbow', abs_levels=None, rel_levels=None, min_abs_level=None,
                  min_rel_level=None, k=2, vinc=2, contours_mask=None, colors_mask=None,
                  vectors_mask=None, color_clim=None, outfile=None, outdir=None, close=False,
                  colorbar_label=None, show=True, contour_color='k', vector_color="k", plot_colorbar=True,
                  max_vector_value_length=5., mas_in_pixel=None, vector_enlarge_factor=1.0,
                  label_size=14, figsize=(20, 5), fig=None, contour_linewidth=1.0, quiver_linewidth=1.0, plot_title=None):
    """
    :param contours: (optional)
        Numpy 2D array (possibly masked) that should be plotted using contours.
    :param colors: (optional)
        Numpy 2D array (possibly masked) that should be plotted using colors.
    :param vectors: (optional)
        Numpy 2D array (possibly masked) that should be plotted using vectors.
    :param vectors_values: (optional)
        Numpy 2D array (possibly masked) that should be used as vector's lengths
        when plotting ``vectors`` array.
    :param cmap: (optional)
        Colormap to use for plotting colors.
        (default: ``gist_rainbow``)
    :param abs_levels: (optional)
        Iterable of absolute levels. If ``None`` then construct levels in other
        way. (default: ``None``)
    :param min_abs_level: (optional)
        Values of minimal absolute level. Used with conjunction of ``k``
        argument for building sequence of absolute levels. If ``None`` then
        construct levels in other way. (default: ``None``)
    :param rel_levels: (optional)
        Iterable of relative levels. If ``None`` then construct levels in other
        way. (default: ``None``)
    :param min_rel_level: (optional)
        Values of minimal relative level. Used with conjunction of ``k``
        argument for building sequence of relative levels. If ``None`` then
        construct levels in other way. (default: ``None``)
    :param k: (optional)
        Factor of incrementation for levels. (default: ``2.0``)
    :param colorbar_label: (optional)
        String to label colorbar. If ``None`` then don't label. (default:
        ``None``)
    :param plot_colorbar: (optional)
        If colors is set then should we plot colorbar? (default: ``True``).
    :param max_vector_value_length: (optional)
        Determines what part of the image is the length of the vector with
        maximum magnitude. E.g. if ``5`` then maximum value of vector quantity
        corresponds to arrow with length equal to 1/5 of the image length.
        (default: ``5``)
    :param mas_in_pixel: (optonal)
        Number of milliarcseconds in one pixel. If ``None`` then plot in pixels.
        (default: ``None``)
    :param vector_enlarge_factor: (optional)
        Additional factor to increase length of vectors representing direction and values of linear polarization.
    """
    matplotlib.rcParams['xtick.labelsize'] = label_size
    matplotlib.rcParams['ytick.labelsize'] = label_size
    matplotlib.rcParams['axes.titlesize'] = label_size
    matplotlib.rcParams['axes.labelsize'] = label_size
    matplotlib.rcParams['font.size'] = label_size
    matplotlib.rcParams['legend.fontsize'] = label_size
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42

    image = None
    if contours is not None:
        image = contours
    elif colors is not None and image is None:
        image = colors
    elif vectors is not None and image is None:
        image = vectors
    else:
        raise Exception("No image to plot")

    x = np.arange(image.shape[0]) - image.shape[0]/2
    y = np.arange(image.shape[1]) - image.shape[1]/2
    if mas_in_pixel is not None:
        x *= mas_in_pixel
        y *= mas_in_pixel

    # Optionally mask arrays
    if contours is not None and contours_mask is not None:
        contours = np.ma.array(contours, mask=contours_mask)
    if colors is not None and colors_mask is not None:
        colors = np.ma.array(colors, mask=colors_mask)
    if vectors is not None and vectors_mask is not None:
        vectors = np.ma.array(vectors, mask=vectors_mask)

    # Actually plotting
    if fig is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    else:
        ax = fig.get_axes()[0]

    if plot_title:
        title = ax.set_title(plot_title, fontsize='large')
    # Plot contours
    if contours is not None:
        if abs_levels is None:
            max_level = np.nanmax(contours)
            if rel_levels is not None:
                abs_levels = [-max_level] + [max_level * i for i in rel_levels]
            else:
                if min_abs_level is not None:
                    n_max = int(math.ceil(math.log(max_level / min_abs_level, k)))
                elif min_rel_level is not None:
                    min_abs_level = min_rel_level * max_level / 100.
                    n_max = int(math.ceil(math.log(max_level / min_abs_level, k)))
                else:
                    raise Exception("Not enough information for levels")
                abs_levels = [-min_abs_level] + [min_abs_level * k ** i for i in
                                                 range(n_max)]
        co = ax.contour(y, x, contours, abs_levels, colors=contour_color, linewidths=contour_linewidth)
    if colors is not None:
        im = ax.imshow(colors, interpolation='none',
                       origin='lower', extent=[y[0], y[-1], x[0], x[-1]],
                       cmap=plt.get_cmap(cmap), clim=color_clim)
    if vectors is not None:
        if vectors_values is not None:
            u = vectors_values * np.cos(vectors)
            v = vectors_values * np.sin(vectors)
            max_vector_value = np.max(np.abs(vectors_values))
            scale = max_vector_value_length*max_vector_value/vector_enlarge_factor
        else:
            u = np.cos(vectors)
            v = np.sin(vectors)
            scale = None

        if vectors_mask is not None:
            u = np.ma.array(u, mask=vectors_mask)
            v = np.ma.array(v, mask=vectors_mask)

        vec = ax.quiver(y[::vinc], x[::vinc], u[::vinc, ::vinc],
                        v[::vinc, ::vinc], angles='uv',
                        units='width', headwidth=0., headlength=0., scale=scale,
                        width=0.001, headaxislength=0., pivot='middle',
                        scale_units='width', color=vector_color, linewidths=quiver_linewidth, edgecolors='k')

    # Set equal aspect
    ax.set_aspect('auto')

    if colors is not None:
        if plot_colorbar:
            from mpl_toolkits.axes_grid1 import make_axes_locatable
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.00)
            cb = fig.colorbar(im, cax=cax)
            if colorbar_label is not None:
                cb.set_label(colorbar_label)

    # Saving output
    if outfile:
        if outdir is None:
            outdir = '.'
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        path = os.path.join(outdir, outfile)
        plt.savefig("{}.png".format(path), bbox_inches='tight', dpi=300)

    if show:
        plt.ioff()
        plt.show()
    if close:
        plt.close()

    return fig


def beta(Gamma):
    """
    Velocity in units of speed of light [c].

    :param Gamma:
        Lorentz factor.
    """
    return np.sqrt(Gamma**2.-1.)/Gamma


def delta(Gamma, theta):
    """
    Doppler factor

    :param Gamma:
        Lorentz factor.
    :param theta:
        LOS angle [rad].
    """
    return 1./(Gamma*(1.-beta(Gamma)*np.cos(theta)))


def theta_plasma(theta_obs, Gamma):
    return np.arctan((np.sin(theta_obs)*np.sqrt(1 - beta(Gamma)**2))/(np.cos(theta_obs) - beta(Gamma)))


# NOTE: My PA is "-"this PA
def PA(phi, theta_los, delta):
    sqrt_part = np.sqrt(np.cos(phi)**2 + np.cos(theta_los)**2*np.sin(phi)**2)
    return np.arctan(np.cos(theta_los)**2*np.sin(phi)*np.tan(delta) / (sqrt_part*(np.sin(theta_los) - np.cos(theta_los)**2*np.cos(phi)*np.tan(delta)/sqrt_part)))


def theta_los_rot(phi, theta_los, delta):
    return np.arccos(np.sin(delta)*np.cos(phi)*np.sin(theta_los) + np.cos(delta)*np.cos(theta_los))


def example_changing_theta_plasma():
    matplotlib.use("TkAgg")
    theta_los = np.deg2rad(0.75)
    delta = np.deg2rad(0.2)
    Gamma = 5

    theta_los = np.deg2rad(5)
    delta = np.deg2rad(1.25)
    Gamma = 10

    phi = np.linspace(0, 2*np.pi, 40)
    pa_single_epoch = -PA(phi, theta_los, delta)
    theta_los_single_epoch = theta_los_rot(phi, theta_los, delta)
    theta_plasma_single_epoch = theta_plasma(theta_los_single_epoch, Gamma)

    fig, axes = plt.subplots(3, 1, sharex=True)
    fig.suptitle(r"$\theta_{\rm LOS} = $" + str(np.round(np.rad2deg(theta_los), 2)) + r"$^{\circ}$, $\delta_{\rm LOS}$ = " +
                 str(np.round(np.rad2deg(delta), 2)) + r"$^{\circ}$")
    axes[0].plot(np.rad2deg(phi), np.rad2deg(pa_single_epoch), lw=2)
    axes[1].plot(np.rad2deg(phi), np.rad2deg(theta_los_single_epoch), lw=2)
    axes[2].plot(np.rad2deg(phi), np.rad2deg(theta_plasma_single_epoch), lw=2)
    axes[2].set_xlabel("Precession azimuth angle, deg")
    axes[2].set_xlim([0, 360])
    axes[0].set_ylabel(r"PA, $^{\circ}$")
    axes[1].set_ylabel(r"$\theta_{\rm LOS}$, $^{\circ}$")
    axes[2].set_ylabel(r"$\theta_{\rm plasma}$, $^{\circ}$")
    plt.show()


def generate_model_images(parallels_run_file, cone_half_angle, LOS_angels_rad, epochs, exec_dir, calculon=False):
    cwd = os.getcwd()
    # Construct params file
    with open(f"{parallels_run_file}", "w+") as fo:
        for epoch, los_angle_rad in zip(epochs, LOS_angels_rad):
            fo.write("{} {} {}".format(los_angle_rad, cone_half_angle, epoch))
            fo.write("\n")

    os.chdir(exec_dir)
    n_jobs = 4
    if calculon:
        n_jobs = 44
    os.system("parallel --files --results generate_model_images_epoch_{3}" + f" --joblog log --jobs {n_jobs} -a {parallels_run_file} -n 1 -m --colsep ' ' \"./bk_transfer\"")
    os.chdir(cwd)


if __name__ == "__main__":

    # NOTE: B and N-field models, as well as \Gamma are specified in main.cpp. Thus, after it you want to ``make``.

    # Work on calculon?
    calculon = True
    # Set working directory according to this:
    # FIXME: You should change this accordingly!
    if calculon:
        # Path to the repo
        base_dir = "/home/ilya/github/wandering-jet"
    else:
        base_dir = "/home/ilya/github/time_machine/bk_transfer"

    # If ``None`` -> using all available epochs. Otherwise, use only first ``n_epochs`` epochs. Useful for debugging
    # purposes, when you want to set e.g. ``n_epochs = 2``.
    n_epochs = None

    # Will be used in folder name containing results. Just to distinguish the results obtained with different models
    # of magnetic field or particle density. E.g. ``toroidal``, ``equipartition_toroidal``, ...
    short_model_description = "test"

    # Precession model #################################################################################################
    # LOS angle of the precession axis
    los_ange_deg_0 = 2.0
    # Amplitude of the precession [deg].
    delta_deg = 0.25
    # Period of the precession [years]
    T_years = 20.
    # Initial phase [0, 2*pi]
    phi_0_rad = 0.0
    ####################################################################################################################

    # Jet cone half-angle [deg]
    cone_half_angle_deg = 0.25

    ####################################################################################################################
    ############################# No need to change anything below this line ###########################################
    ####################################################################################################################

    if calculon:
        n_jobs = 44
    else:
        n_jobs = 4

    stokes = ("I", "Q", "U")

    # Directory with pics and UVFITS-files.
    data_dir = "/home/ilya/Downloads/1641+399"
    files = glob.glob(os.path.join(data_dir, "*.pcn.png"))
    times, epochs = get_epochs_from_image_list(files)

    los_angle_rad_0 = np.deg2rad(los_ange_deg_0)
    delta = np.deg2rad(delta_deg)

    cone_half_angle = np.deg2rad(cone_half_angle_deg)

    # Directory to save results
    save_dir = "{}/results/3C454.3/{}_los_{}_delta_{}_T_{}_phi0_{}".format(base_dir, short_model_description, los_ange_deg_0,
                                                                           delta_deg, T_years, phi_0_rad)
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    freq_ghz = 15.4
    # Some template UVFITS with full polarization. Its uv-coverage and noise will be used while creating fake data
    # template_uvfits = "/home/ilya/github/time_machine/bk_transfer/uvfits/1458+718.u.2006_09_06.uvf"
    # Multiplicative factor for noise added to model visibilities.
    noise_scale_factor = 1.0
    # Used in CLEAN
    mapsize = (512, 0.1)

    # 1641+399 stack beam from Pushkarev+2023
    common_beam = (0.73, 0.73, 0)

    npixels_beam = np.pi*common_beam[0]*common_beam[1]/(4*np.log(2)*mapsize[1]**2)
    print("#pixels in beam = {}".format(int(npixels_beam)))

    jetpol_run_directory = "{}/Release".format(base_dir)

    # C++ code run parameters
    # M87
    # z = 0.00436
    # 1641+399
    z = 0.59
    n_along = 500
    n_across = 300
    lg_pixel_size_mas_min = np.log10(0.01)
    lg_pixel_size_mas_max = np.log10(0.1)
    resolutions = np.logspace(lg_pixel_size_mas_min, lg_pixel_size_mas_max, n_along)
    print("Model jet extends up to {:.1f} mas!".format(np.sum(resolutions)))

    # Plot only jet emission and do not plot counter-jet?
    jet_only = True
    path_to_script = "{}/scripts/script_clean_rms".format(base_dir)
    parallels_run_file = "{}/parallels_run.txt".format(base_dir)

    images_i = list()
    images_q = list()
    images_u = list()
    images_pang = list()

    # Find phases, LOS-angles and PA-angles for given period and initial phase
    phis = list()
    LOS_angels_rad = list()
    PAs_rad = list()
    epochs_to_calculate = list()
    times_to_calculate = list()
    for i_epoch, (time, epoch) in enumerate(zip(times, epochs)):
        if n_epochs is not None:
            if i_epoch == n_epochs:
                break
        # phi = np.random.uniform(0, 2*np.pi, 1)[0]
        phi = phi_0_rad + (time % T_years)/T_years*2*np.pi
        if phi > 2*np.pi:
            phi -= 2*np.pi
        phis.append(phi)
        print("Phase = {:.2f} rad".format(phi))
        pa_single_epoch = -PA(phi, los_angle_rad_0, delta)
        PAs_rad.append(pa_single_epoch)
        los_angle_rad = theta_los_rot(phi, los_angle_rad_0, delta)
        LOS_angels_rad.append(los_angle_rad)
        epochs_to_calculate.append(epoch)
        times_to_calculate.append(time)

    # Calculate only ``n_epochs`` epochs.
    epochs = epochs_to_calculate
    times = times_to_calculate
    # Generate images in parallel
    generate_model_images(parallels_run_file, cone_half_angle, LOS_angels_rad, epochs, jetpol_run_directory,
                          calculon=calculon)

    rot_angles_deg = list()
    # Make true pics, generate and CLEAN single epoch uv-data
    for i_epoch, (time, epoch) in enumerate(zip(times, epochs)):
        print("Time = {:.2f} yrs, epoch = {}".format(time, epoch))
        template_uvfits = os.path.join(data_dir, "1641+399.u.{}.uvf".format(epoch))

        los_angle_rad = LOS_angels_rad[i_epoch]
        pa_single_epoch = PAs_rad[i_epoch]

        # Jet to the West as in 1641+399
        rot_angle_deg = -90 + np.rad2deg(pa_single_epoch)
        rot_angles_deg.append(rot_angle_deg)

        print("LOS angle (deg) = {:.2f}\n".format(np.rad2deg(los_angle_rad)))
        print("PA angle (deg) = {:.2f}\n".format(rot_angle_deg))

        # Plot true image of polarization
        imagei = np.loadtxt(os.path.join(jetpol_run_directory, "{}/jet_image_{}_{}_{}.txt".format(jetpol_run_directory, "i", freq_ghz, epoch)))
        imageq = np.loadtxt(os.path.join(jetpol_run_directory, "{}/jet_image_{}_{}_{}.txt".format(jetpol_run_directory, "q", freq_ghz, epoch)))
        imageu = np.loadtxt(os.path.join(jetpol_run_directory, "{}/jet_image_{}_{}_{}.txt".format(jetpol_run_directory, "u", freq_ghz, epoch)))
        mask = imagei == 0
        imagei[mask] = np.nan
        imageq[mask] = np.nan
        imageu[mask] = np.nan
        imagep = np.hypot(imageq, imageu)
        imagepang = 0.5*np.arctan2(imageu, imageq)

        min_abs_lev = 0.001*np.max(imagei)

        fig = plot_function(contours=imagei, colors=imagep, vectors=imagepang,
                            vectors_values=None, min_rel_level=0.001,
                            vinc=10, contour_color="gray", vector_color="k", cmap="gist_rainbow",
                            vector_enlarge_factor=8, colorbar_label="PPOL, Jy/pixel", contour_linewidth=1.0,
                            quiver_linewidth=0.01, plot_title="LOS = {:.2f} deg, rot = {:.2f} deg".format(np.rad2deg(los_angle_rad), rot_angle_deg))
        fig = plot_function(contours=imagei, abs_levels=[0.01*np.max(imagei)], fig=fig, show=False, close=True)
        fig.savefig(os.path.join(save_dir, "true_pol_{}.png".format(epoch)), dpi=300, bbox_inches="tight")
        plt.close()




    epochs_string = " ".join([epoch for epoch in epochs])
    rot_angles_string = " ".join([str(rot) for rot in rot_angles_deg])
    os.system("parallel --link -k --results create_artificial_epoch_{1} --jobs %d \"python %s/create_artificial.py --save_dir %s --data_dir %s  --exec_dir %s"
              " --epoch {1} --rot_angle_deg {2} \" ::: %s ::: %s" % (n_jobs, base_dir, save_dir, data_dir,
                                                                     jetpol_run_directory, epochs_string,
                                                                     rot_angles_string))

        # # Substitute real uv-data with model value, rotate jet and polarization ########################################
        # uvdata = UVData(template_uvfits)
        # noise = uvdata.noise(average_freq=False, use_V=False)
        # # If one needs to decrease the noise this is the way to do it
        # for baseline, baseline_noise_std in noise.items():
        #     noise.update({baseline: noise_scale_factor*baseline_noise_std})
        #
        # stokes = ("I", "Q", "U", "V")
        # jms = [JetImage(z=z, n_along=n_along, n_across=n_across,
        #                 lg_pixel_size_mas_min=lg_pixel_size_mas_min, lg_pixel_size_mas_max=lg_pixel_size_mas_max,
        #                 jet_side=True, rot=np.deg2rad(rot_angle_deg)) for _ in stokes]
        # # cjms = [JetImage(z=z, n_along=n_along, n_across=n_across,
        # #                  lg_pixel_size_mas_min=lg_pixel_size_mas_min, lg_pixel_size_mas_max=lg_pixel_size_mas_max,
        # #                  jet_side=False) for _ in stokes]
        # for i, stk in enumerate(stokes):
        #     jms[i].load_image_stokes(stk, "{}/jet_image_{}_{}_{}.txt".format(jetpol_run_directory, stk.lower(), freq_ghz, epoch), scale=1.0)
        #     # cjms[i].load_image_stokes(stk, "../{}/cjet_image_{}_{}.txt".format(jetpol_run_directory, stk.lower(), freq_ghz), scale=1.0)
        #
        # # List of models (for J & CJ) for all stokes
        # # js = [TwinJetImage(jms[i], cjms[i]) for i in range(len(stokes))]
        #
        # uvdata.zero_data()
        # if jet_only:
        #     uvdata.substitute(jms)
        # else:
        #     # uvdata.substitute(js)
        #     pass
        # # Rotate EVPA also
        # uvdata.rotate_evpa(np.deg2rad(rot_angle_deg))
        # uvdata.noise_add(noise)
        # if epoch in ("2014_06_05",):
        #     downscale_by_freq = True
        # else:
        #     downscale_by_freq = False
        # uvdata.save(os.path.join(save_dir, "artificial_{}.uvf".format(epoch)), rewrite=True, downscale_by_freq=downscale_by_freq)









        # # CLEAN synthetic UV-data
        # for stk in stokes:
        #     outfname = "model_cc_{}_{}.fits".format(epoch, stk.lower())
        #     if os.path.exists(os.path.join(save_dir, outfname)):
        #         os.unlink(os.path.join(save_dir, outfname))
        #     clean_difmap(fname="template.uvf", path=save_dir,
        #                  outfname=outfname, outpath=save_dir, stokes=stk.lower(),
        #                  mapsize_clean=mapsize, path_to_script=path_to_script,
        #                  show_difmap_output=False,
        #                  # text_box=text_boxes[freq],
        #                  dfm_model=os.path.join(save_dir, "model_cc_{}.mdl".format(stk.lower())),
        #                  beam_restore=common_beam)


    # CLEAN synthetic UV-data in parallel
    fnames = " ".join(["artificial_{}.uvf".format(epoch) for epoch in epochs])
    os.system(f"parallel -k --jobs {n_jobs} python {base_dir}/clean_uvfits.py --mapsize_clean {mapsize[0]} {mapsize[1]} "
              f"--save_dir \"{save_dir}\" --path_to_script \"{path_to_script}\"  --beam_restore {common_beam[0]} {common_beam[1]} {common_beam[2]}  "
              f"--fname ::: {fnames}")

    # Plot pictures and make stack
    for i_epoch, (time, epoch) in enumerate(zip(times, epochs)):
        ccimages = {stk: create_clean_image_from_fits_file(os.path.join(save_dir, "artificial_{}_{}.fits".format(epoch, stk.lower())))
                    for stk in stokes}
        ipol = ccimages["I"].image
        beam = ccimages["I"].beam
        # Number of pixels in beam

        std = find_image_std(ipol, beam_npixels=npixels_beam)
        print("IPOL image std = {} mJy/beam".format(1000*std))
        blc, trc = find_bbox(ipol, level=4*std, min_maxintensity_mjyperbeam=100*std,
                             min_area_pix=20*npixels_beam, delta=10)
        if blc[0] == 0: blc = (blc[0]+1, blc[1])
        if blc[1] == 0: blc = (blc[0], blc[1]+1)
        if trc[0] == ipol.shape: trc = (trc[0]-1, trc[1])
        if trc[1] == ipol.shape: trc = (trc[0], trc[1]-1)
        masks_dict, ppol_quantile = pol_mask({stk: ccimages[stk].image for stk in stokes}, npixels_beam, n_sigma=4,
                                             return_quantile=True)
        ppol = np.hypot(ccimages["Q"].image, ccimages["U"].image)
        ppol = correct_ppol_bias(ipol, ppol, ccimages["Q"].image, ccimages["U"].image, npixels_beam)
        pang = 0.5*np.arctan2(ccimages["U"].image, ccimages["Q"].image)
        fpol = ppol/ipol

        images_i.append(ipol)
        images_q.append(ccimages["Q"].image)
        images_u.append(ccimages["U"].image)
        images_pang.append(np.ma.array(pang, mask=masks_dict["P"]))

        # Make a single epoch map
        # PPOL contours
        fig = iplot(ppol, x=ccimages["I"].x, y=ccimages["I"].y,
                    min_abs_level=ppol_quantile, blc=blc, trc=trc,
                    close=False, contour_color='black',
                    plot_colorbar=False)
        # Add single IPOL contour and vectors of the PANG
        fig = iplot(contours=ipol, vectors=pang,
                    x=ccimages["I"].x, y=ccimages["I"].y, vinc=4, contour_linewidth=0.25,
                    vectors_mask=masks_dict["P"], abs_levels=[3*std], blc=blc, trc=trc,
                    beam=common_beam, close=True, show_beam=True, show=False,
                    contour_color='gray', fig=fig, vector_color="black", plot_colorbar=False)
        axes = fig.get_axes()[0]
        axes.invert_xaxis()
        fig.savefig(os.path.join(save_dir, f"single_epoch_p_epoch_{epoch}.png"), dpi=600, bbox_inches="tight")
        plt.close()

    np.savetxt(os.path.join(save_dir, "phis.txt"), phis)

    stack_i = np.mean(images_i, axis=0)
    stack_q = np.mean(images_q, axis=0)
    stack_u = np.mean(images_u, axis=0)
    stack_p = np.hypot(stack_q, stack_u)
    stack_pang = 0.5*np.arctan2(stack_u, stack_q)
    stack_pang_std = stat_of_masked(images_pang, stat="scipy_circstd", n_epochs_not_masked_min=3)
    stack_p = correct_ppol_bias(stack_i, stack_p, stack_q, stack_u, npixels_beam)
    stack_fpol = stack_p/stack_i
    masks_dict, ppol_quantile = pol_mask({"I": stack_i, "Q": stack_q, "U": stack_u}, npixels_beam, n_sigma=4,
                                         return_quantile=True)
    std = find_image_std(stack_i, beam_npixels=npixels_beam)
    blc, trc = find_bbox(stack_i, level=4*std, min_maxintensity_mjyperbeam=100*std,
                         min_area_pix=20*npixels_beam, delta=10)

    # IPOL contours
    fig = iplot(stack_i, x=ccimages["I"].x, y=ccimages["I"].y,
                min_abs_level=4*std, blc=blc, trc=trc, beam=beam, close=True, show_beam=True, show=False,
                contour_color='gray', contour_linewidth=0.25)
    fig.savefig(os.path.join(save_dir, "observed_stack_i.png"), dpi=600, bbox_inches="tight")

    # PPOL contours
    fig = iplot(stack_p, x=ccimages["I"].x, y=ccimages["I"].y,
                min_abs_level=ppol_quantile, blc=blc, trc=trc,
                close=False, contour_color='black',
                plot_colorbar=False)
    # Add single IPOL contour and vectors of the PANG
    fig = iplot(contours=stack_i, vectors=stack_pang,
                x=ccimages["I"].x, y=ccimages["I"].y, vinc=4, contour_linewidth=0.25,
                vectors_mask=masks_dict["P"], abs_levels=[3*std], blc=blc, trc=trc,
                beam=common_beam, close=True, show_beam=True, show=False,
                contour_color='gray', fig=fig, vector_color="black", plot_colorbar=False)
    axes = fig.get_axes()[0]
    axes.invert_xaxis()
    fig.savefig(os.path.join(save_dir, "observed_stack_p.png"), dpi=600, bbox_inches="tight")

    fig = iplot(stack_i, stack_fpol, x=ccimages["I"].x, y=ccimages["I"].y,
                min_abs_level=4*std, colors_mask=masks_dict["P"], color_clim=[0, 0.7], blc=blc, trc=trc,
                beam=common_beam, close=True, colorbar_label="m", show_beam=True, show=False,
                cmap='gnuplot', contour_color='black', plot_colorbar=True,
                contour_linewidth=0.25)
    fig.savefig(os.path.join(save_dir, "observed_stack_fpol.png"), dpi=600, bbox_inches="tight")

    fig = iplot(stack_i, np.rad2deg(stack_pang_std), x=ccimages["I"].x, y=ccimages["I"].y,
                min_abs_level=4*std, colors_mask=masks_dict["P"], color_clim=None, blc=blc, trc=trc,
                beam=common_beam, close=True, colorbar_label=r"$\sigma_{\rm EVPA},$ $^{\circ}$", show_beam=True, show=False,
                cmap='gnuplot', contour_color='black', plot_colorbar=True,
                contour_linewidth=0.25)
    fig.savefig(os.path.join(save_dir, "observed_stack_pang_std.png"), dpi=600, bbox_inches="tight")

    # Save stack
    stack_dict = {"stack_i": stack_i, "stack_q": stack_q, "stack_u": stack_u, "stack_p": stack_p,
                  "stack_fpol": stack_fpol, "stack_pang": stack_pang, "stack_pang_std": stack_pang_std}
    for stack_type, stack_array in stack_dict.items():
        np.savetxt(os.path.join(save_dir, "{}.txt".format(stack_type)), stack_array)