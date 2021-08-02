import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
from jet_image import JetImage, TwinJetImage, convert_difmap_model_file_to_CCFITS


# This script gets txt-images produced by C++ code and creates Difmap-format model files and FITS-files at several
# frequencies. Common beam and mapsize must be specified.

# Directory to save files
save_dir = "/home/ilya/github/bk_transfer/results"
# Observed frequencies (GHz) of simulations. Must coincide with that in C++ code
# freqs_obs_ghz = [8.1, 12.1, 15.4]
freqs_obs_ghz = [15.4]
# Will be used for creating FITS images at each frequency
common_mapsize = (1024, 0.1)
# Common beam in (mas, mas, deg). Can be arbitrary small/big
common_beam = (1, 1, 0)
# common_beam = None
# Where C++ result txt files?
jetpol_run_directory = "Release"
# C++ code run parameters
z = 0.1
n_along = 600
n_across = 400
lg_pixel_size_mas_min = np.log10(0.01)
lg_pixel_size_mas_max = np.log10(0.1)

##############################################
# No need to change anything below this line #
##############################################

# Plot only jet emission and do not plot counter-jet?
jet_only = True
# For some reason Difmap needs UVFITS file (to make invert to build a map). But for creating FITS image anyone could be
# used
template_uvfits_dict = {15.4: "1458+718.u.2006_09_06.uvf", 12.1: "1458+718.j.2006_09_06.uvf",
                        8.4: "1458+718.y.2006_09_06.uvf", 8.1: "1458+718.x.2006_09_06.uvf"}
template_uvfits_dir = "/home/ilya/github/bk_transfer/uvfits"

common_mapsize_x2 = (int(common_mapsize[0]*2), common_mapsize[1])

for freq in freqs_obs_ghz:
    stokes = ("I", "Q", "U")
    jms = [JetImage(z=z, n_along=n_along, n_across=n_across,
                    lg_pixel_size_mas_min=lg_pixel_size_mas_min, lg_pixel_size_mas_max=lg_pixel_size_mas_max,
                    jet_side=True) for _ in stokes]
    # cjms = [JetImage(z=z, n_along=n_along, n_across=n_across,
    #                  lg_pixel_size_mas_min=lg_pixel_size_mas_min, lg_pixel_size_mas_max=lg_pixel_size_mas_max,
    #                  jet_side=False) for _ in stokes]
    for i, stk in enumerate(stokes):
        jms[i].load_image_stokes(stk, "{}/jet_image_{}_{}.txt".format(jetpol_run_directory, stk.lower(), freq))
        # cjms[i].load_image_stokes(stk, "{}/cjet_image_{}_{}.txt".format(jetpol_run_directory, stk.lower(), freq))

    # List of models (for J & CJ) for all stokes
    # js = [TwinJetImage(jms[i], cjms[i]) for i in range(len(stokes))]

    template_uvfits = os.path.join(template_uvfits_dir, template_uvfits_dict[freq])

    if jet_only:
        for stk, jm in zip(stokes, jms):
            print("======================================================================================")
            print("      Saving model image (Stokes {}, frequency {} GHz) to difmap format...".format(stk, freq))
            print("======================================================================================")
            jm.save_image_to_difmap_format(os.path.join(save_dir, "model_dfm_{}_{}.mdl".format(stk.lower(), freq)))
            print("======================================================================================")
            print("     Convolving with common beam ({} mas, {} mas, {} deg) and saving to FITS file...".format(*common_beam))
            print("======================================================================================")
            convert_difmap_model_file_to_CCFITS(os.path.join(save_dir, "model_dfm_{}_{}.mdl".format(stk.lower(), freq)),
                                                stk, common_mapsize_x2, common_beam, template_uvfits,
                                                os.path.join(save_dir, "convolved_{}_{}.fits".format(stk.lower(), freq)))
    else:
        raise NotImplementedError
