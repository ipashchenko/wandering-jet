import os
import sys
import argparse
import numpy as np
from jet_image import JetImage, TwinJetImage
sys.path.insert(0, 've/vlbi_errors')
from uv_data import UVData


if __name__ == "__main__":

    # Substitute real uv-data with model value, rotate jet and polarization ########################################

    CLI = argparse.ArgumentParser()
    CLI.add_argument("--epoch", type=str)
    CLI.add_argument("--rot_angle_deg", type=float)
    CLI.add_argument("--noise_scale", type=float, default=1.0)
    CLI.add_argument("--jet_only", type=bool, default=True)
    CLI.add_argument("--lg_pixel_size_mas_min", type=float, default=-2.)
    CLI.add_argument("--lg_pixel_size_mas_max", type=float, default=-1.)
    CLI.add_argument("--redshift", type=float, default=0.59)
    CLI.add_argument("--nalong", type=int, default=500)
    CLI.add_argument("--nacross", type=int, default=300)
    CLI.add_argument("--save_dir", type=str)
    CLI.add_argument("--data_dir", type=str)
    CLI.add_argument("--exec_dir", type=str)

    args = CLI.parse_args()

    epoch = args.epoch
    noise_scale_factor = args.noise_scale
    save_dir = args.save_dir
    data_dir = args.data_dir
    z = args.redshift
    n_along = args.nalong
    n_across = args.nacross
    lg_pixel_size_mas_min = args.lg_pixel_size_mas_min
    lg_pixel_size_mas_max = args.lg_pixel_size_mas_max
    rot_angle_deg = args.rot_angle_deg
    exec_dir = args.exec_dir
    jet_only = args.jet_only

    print(epoch)
    print(noise_scale_factor)
    print(save_dir)
    print(data_dir)
    print(z)
    print(n_along)
    print(n_across)
    print(lg_pixel_size_mas_min)
    print(lg_pixel_size_mas_max)
    print(rot_angle_deg)
    print(jet_only)

    template_uvfits = os.path.join(data_dir, "1641+399.u.{}.uvf".format(epoch))
    print("Using template UVFITS: ", template_uvfits)

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
        jms[i].load_image_stokes(stk, "{}/jet_image_{}_{}_{}.txt".format(exec_dir, stk.lower(), 15.4, epoch), scale=1.0)
        # cjms[i].load_image_stokes(stk, "../{}/cjet_image_{}_{}.txt".format(jetpol_run_directory, stk.lower(), freq_ghz), scale=1.0)

    # List of models (for J & CJ) for all stokes
    # js = [TwinJetImage(jms[i], cjms[i]) for i in range(len(stokes))]

    uvdata.zero_data()
    if jet_only:
        uvdata.substitute(jms)
    else:
        # uvdata.substitute(js)
        pass
    # Rotate EVPA also
    uvdata.rotate_evpa(np.deg2rad(rot_angle_deg))
    uvdata.noise_add(noise)
    if epoch in ("2014_06_05",):
        downscale_by_freq = True
    else:
        downscale_by_freq = False
    uvdata.save(os.path.join(save_dir, "artificial_{}.uvf".format(epoch)), rewrite=True, downscale_by_freq=downscale_by_freq)


