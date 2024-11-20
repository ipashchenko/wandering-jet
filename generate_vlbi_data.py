import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import sys
from jet_image import JetImage, TwinJetImage
from vlbi_utils import find_image_std, find_bbox, pol_mask, correct_ppol_bias
sys.path.insert(0, '../ve/vlbi_errors')
from uv_data import UVData
from spydiff import clean_difmap
from from_fits import create_clean_image_from_fits_file
from image import plot as iplot

# -107 for M87
rot_angle_deg = 0.0
freq_ghz = 15.4
# Directory to save files
save_dir = "/home/ilya/github/time_machine/bk_transfer/results"
# Some template UVFITS with full polarization. Its uv-coverage and noise will be used while creating fake data
template_uvfits = "/home/ilya/github/time_machine/bk_transfer/uvfits/1458+718.u.2006_09_06.uvf"
# Multiplicative factor for noise added to model visibilities.
noise_scale_factor = 0.1
# Used in CLEAN
mapsize = (512, 0.1)
jetpol_run_directory = "/home/ilya/github/time_machine/bk_transfer/Release"

# C++ code run parameters
# z = 0.00436
z = 0.1
n_along = 1000
n_across = 400
lg_pixel_size_mas_min = np.log10(0.01)
lg_pixel_size_mas_max = np.log10(0.1)
resolutions = np.logspace(lg_pixel_size_mas_min, lg_pixel_size_mas_max, n_along)
print("Model jet extends up to {:.1f} mas!".format(np.sum(resolutions)))
##############################################
# No need to change anything below this line #
##############################################

# Plot only jet emission and do not plot counter-jet?
jet_only = True
path_to_script = "/home/ilya/github/time_machine/bk_transfer/scripts/script_clean_rms"
uvdata = UVData(template_uvfits)
noise = uvdata.noise(average_freq=False, use_V=False)
# If one needs to decrease the noise this is the way to do it
for baseline, baseline_noise_std in noise.items():
    noise.update({baseline: noise_scale_factor*baseline_noise_std})

stokes = ("I", "Q", "U", "V")
jms = [JetImage(z=z, n_along=n_along, n_across=n_across,
                lg_pixel_size_mas_min=lg_pixel_size_mas_min, lg_pixel_size_mas_max=lg_pixel_size_mas_max,
                jet_side=True, rot=np.deg2rad(rot_angle_deg)) for _ in stokes]
# cjms = [JetImage(z=z, n_along=n_along, n_across=n_across,
#                  lg_pixel_size_mas_min=lg_pixel_size_mas_min, lg_pixel_size_mas_max=lg_pixel_size_mas_max,
#                  jet_side=False) for _ in stokes]
for i, stk in enumerate(stokes):
    jms[i].load_image_stokes(stk, "{}/jet_image_{}_{}.txt".format(jetpol_run_directory, stk.lower(), freq_ghz), scale=1.0)
    # cjms[i].load_image_stokes(stk, "../{}/cjet_image_{}_{}.txt".format(jetpol_run_directory, stk.lower(), freq_ghz), scale=1.0)

# List of models (for J & CJ) for all stokes
# js = [TwinJetImage(jms[i], cjms[i]) for i in range(len(stokes))]

uvdata.zero_data()
if jet_only:
    uvdata.substitute(jms)
else:
    # uvdata.substitute(js)
    pass
# Optionally
uvdata.rotate_evpa(np.deg2rad(rot_angle_deg))
uvdata.noise_add(noise)
uvdata.save(os.path.join(save_dir, "template.uvf"), rewrite=True)

# CLEAN synthetic UV-data
for stk in stokes:
    outfname = "model_cc_{}.fits".format(stk.lower())
    if os.path.exists(os.path.join(save_dir, outfname)):
        os.unlink(os.path.join(save_dir, outfname))
    clean_difmap(fname="template.uvf", path=save_dir,
                 outfname=outfname, outpath=save_dir, stokes=stk.lower(),
                 mapsize_clean=mapsize, path_to_script=path_to_script,
                 show_difmap_output=True,
                 # text_box=text_boxes[freq],
                 dfm_model=os.path.join(save_dir, "model_cc_{}.mdl".format(stk.lower())))

ccimages = {stk: create_clean_image_from_fits_file(os.path.join(save_dir, "model_cc_{}.fits".format(stk.lower())))
            for stk in stokes}
ipol = ccimages["I"].image
beam = ccimages["I"].beam
# Number of pixels in beam
npixels_beam = np.pi*beam[0]*beam[1]/(4*np.log(2)*mapsize[1]**2)

std = find_image_std(ipol, beam_npixels=npixels_beam)
print("IPOL image std = {} mJy/beam".format(1000*std))
blc, trc = find_bbox(ipol, level=4*std, min_maxintensity_mjyperbeam=10*std,
                     min_area_pix=4*npixels_beam, delta=10)
if blc[0] == 0: blc = (blc[0]+1, blc[1])
if blc[1] == 0: blc = (blc[0], blc[1]+1)
if trc[0] == ipol.shape: trc = (trc[0]-1, trc[1])
if trc[1] == ipol.shape: trc = (trc[0], trc[1]-1)
masks_dict, ppol_quantile = pol_mask({stk: ccimages[stk].image for stk in
                                     ("I", "Q", "U")}, npixels_beam, n_sigma=4,
                                     return_quantile=True)
ppol = np.hypot(ccimages["Q"].image, ccimages["U"].image)
ppol = correct_ppol_bias(ipol, ppol, ccimages["Q"].image, ccimages["U"].image, npixels_beam)
pang = 0.5*np.arctan2(ccimages["U"].image, ccimages["Q"].image)
fpol = ppol/ipol

# IPOL contours
fig = iplot(ipol, x=ccimages["I"].x, y=ccimages["I"].y,
            min_abs_level=4*std, blc=blc, trc=trc, beam=beam, close=True, show_beam=True, show=False,
            contour_color='gray', contour_linewidth=0.25)
fig.savefig(os.path.join(save_dir, "observed_ipol.png"), dpi=600, bbox_inches="tight")

# PPOL contours
fig = iplot(ppol, x=ccimages["I"].x, y=ccimages["I"].y,
            min_abs_level=ppol_quantile, blc=blc, trc=trc,
            close=False, contour_color='black',
            plot_colorbar=False)
# Add single IPOL contour and vectors of the PANG
fig = iplot(contours=ipol, vectors=pang,
            x=ccimages["I"].x, y=ccimages["I"].y, vinc=4, contour_linewidth=0.25,
            vectors_mask=masks_dict["P"], abs_levels=[3*std], blc=blc, trc=trc,
            beam=beam, close=True, show_beam=True, show=False,
            contour_color='gray', fig=fig, vector_color="black", plot_colorbar=False)
axes = fig.get_axes()[0]
axes.invert_xaxis()
fig.savefig(os.path.join(save_dir, "observed_ppol.png"), dpi=600, bbox_inches="tight")

fig = iplot(ipol, fpol, x=ccimages["I"].x, y=ccimages["I"].y,
            min_abs_level=4*std, colors_mask=masks_dict["P"], color_clim=[0, 0.7], blc=blc, trc=trc,
            beam=beam, close=True, colorbar_label="m", show_beam=True, show=False,
            cmap='gnuplot', contour_color='black', plot_colorbar=True,
            contour_linewidth=0.25)
fig.savefig(os.path.join(save_dir, "observed_fpol.png"), dpi=600, bbox_inches="tight")
