import os
import sys
import argparse
sys.path.insert(0, '/home/ilya/github/ve/vlbi_errors')
from spydiff import CLEAN_difmap, clean_difmap



os.system(f"parallel -k --jobs {n_jobs} python {script_dir}/clean_uvfits.py --beam_fractions \"{beam_fracs}\" --mapsize_clean \"{mapsize} {mapsizes_dict[freq_ghz][1]}\" --save_dir \"{save_dir}\" --path_to_script \"{path_to_script}\"  --nw_beam_size \"{nw_beam_size}\" --use_elliptical \"{use_elliptical}\" --fname ::: {fnames}")

if __name__ == "__main__":

    CLI = argparse.ArgumentParser()
    CLI.add_argument("--fname",
                     type=str)
    # CLI.add_argument("--beam_fractions",  # name on the CLI - drop the `--` for positional/required parameters
    #                  nargs="*",  # 0 or more values expected => creates a list
    #                  type=float,
    #                  default=[1.0],  # default if nothing is provided
    #                  )
    CLI.add_argument("--mapsize_clean",  # name on the CLI - drop the `--` for positional/required parameters
                     nargs="*",  # 0 or more values expected => creates a list
                     type=float,
                     default=[512, 0.1],  # default if nothing is provided
                     )
    CLI.add_argument("--beam_restore",  # name on the CLI - drop the `--` for positional/required parameters
                     type=float
                     # default=[1.353, 1.557, -65.65],  # default if nothing is provided
                     # default=[4.561, 5.734, -51.67],  # default if nothing is provided
                     )
    CLI.add_argument("--save_dir",
                     type=str)
    CLI.add_argument("--path_to_script",
                     type=str)
    args = CLI.parse_args()

    uvfits_file = args.fname
    save_dir = args.save_dir
    beam_restore = args.beam_restore
    path_to_script = args.path_to_script
    mapsize_clean = args.mapsize_clean

    stokes = ("I", "Q", "U")
    epoch = None
    for stk in stokes:
        outfname = "model_cc_{}_{}.fits".format(epoch, stk.lower())
        if os.path.exists(os.path.join(save_dir, outfname)):
            os.unlink(os.path.join(save_dir, outfname))
        clean_difmap(fname=uvfits_file, path=save_dir,
                     outfname=outfname, outpath=save_dir, stokes=stk.lower(),
                     mapsize_clean=mapsize_clean, path_to_script=path_to_script,
                     show_difmap_output=False,
                     beam_restore=beam_restore)
        # CLEAN_difmap(uvfits_file, stk.lower(), mapsize, outfname, restore_beam=beam_restore,
        #              boxfile=None, working_dir=None, uvrange=None,
        #              box_clean_nw_niter=1000, clean_gain=0.03, dynam_su=20, dynam_u=6, deep_factor=1.0,
        #              remove_difmap_logs=True, save_noresid=None, save_resid_only=None, save_dfm=None,
        #              noise_to_use="F", shift=None)
