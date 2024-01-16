# carp
CCD Astronomy Reduction Pipeline
This is the CCD Astro Reduction Pipeline (CARP).
Authors:
    Midge Hartshorn
    Sophie Booth

BASE FUNCTIONS:
    collect_files(datadir, cal_status)
    median_combine(filelist)
    bias_subtract(filename, path_to_master_bias, output_path)
    scale_master_dark(path_to_master_dark, desired_exptime)
    dark_subtract(filename, path_to_master_dark, outpath)
    norm_combine_flats(filelist)
    flatfield(filename, path_to_master_flat, outpath)
AUTOMATED FUNCTIONS:
    create_master_bias(biaspath)
    create_master_dark(darkpath, biaspath, exposuretime)
    create_master_flat(flatpath, darkpath, biaspath, targ, filter, exposuretime)
    reduce_images(datapath, biaspath, darkpath, flatpath, objectname, filtername, exposuretime)
    align_and_stack()
