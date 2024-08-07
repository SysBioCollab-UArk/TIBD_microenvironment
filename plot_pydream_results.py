from run_tibd_pydream import *  # import from the file you ran PyDREAM from
import glob
import os

# Create plots using output files already generated by PyDREAM using the calibrator.create_figures() function

# get the existing PyDREAM output files
path = os.getcwd()  # path to where PyDREAM generated files are (default is current working diretory)
logps_files = glob.glob(os.path.join(path, 'dreamzs*logps*'))  # need to pass ALL 'logps' files
samples_files = glob.glob(os.path.join(path, 'dreamzs*params*'))  # need to pass ALL 'params' files

# create the ParameterCalibration object (copy from the file you ran PyDREAM from)
calibrator = ParameterCalibration(model,
                                  exp_data_file,
                                  [tumor_injection] * 2 + [tumor_bisphos_injection, bisphos_injection],
                                  priors=custom_priors,
                                  no_sample=no_sample,
                                  param_expts_map=param_expts_map)

# call the 'create_figures' function
calibrator.create_figures(logps_files, samples_files, show_plots=True, plot_tc_args={'separate_plots': True})
