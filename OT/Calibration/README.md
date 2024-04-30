Codes to process .txt files from applied waves and calibrate the trap using the results. Relies on other codes/functions not found here (these are just the ones I made significant changes to).

##### _Step1_get_TF_PSx_PSy_all_files.m_
- Modified version of loic_codes/TF/Step1_TF_Loic_2018b_read_all_using_function.m
- Calculates transfer function & power spectrum from .txt files from multiharmonic OT measurements
- Biggest difference is that it calculates the power spectrum in both x & y (parallel & perpendicular to trap mvmt)
- A few other changes you can check out in the code itself
- Make sure to set the Fexc to match the wave you were using! Unless you're me, it's probably not the ones from "Ora's fun experiments"
- RELIES ON LOTS OF OTHER FUNCTIONS NOT FOUND IN THIS FOLDER!
- **requires _oras_codes/Optical Trap Calibration/Step1_get_TF_PSx_PSy_FUNCTION.m_, _loic_codes/TF/power_spectrum_Adam_Loic_FUNCTION.m_, _loic_codes/TF/loic_detect_changepoint_QPD_FUNCTION.m_, _loic_codes/TF/NFFT_generator_loic_v4_constant_windows.m_**

##### _Step1_get_PSy_all_files.m_
- Modified version of Step1_TF_Loic_2018b_read_all_using_function.m
- Only calculates the power spectrum perpendicular to trap mvmt from .txt files from multiharmonic OT measurements
- Useful if you've already used Step1_TF_Loic_2018b_read_all_using_function.m but now also want the yPS
- RELIES ON OTHER FUNCTIONS NOT FOUND IN THIS FOLDER!
- **requires _loic_codes/TF/power_spectrum_Adam_Loic_FUNCTION.m_**

##### _OTC_sqwv_beta_estimate.m_
- Modified version of loic_codes/adams-method/Fitting_adam_method_r10.m to carry out optical trap calibration (OTC - geddit?). 
- Takes .mats generated from loic_codes/TF/Step1_TF_Loic_2018b_read_all_using_function.m or oras_codes/Optical Trap Calibration/Step1_get_TF_PSx_PSy_all_files.m and fits magnitude, phase and power spectrum to get relevant parameters
  - namely Beta, GammaR, Alpha, k_trap, k_cyt0, k_cyt1
- Requires a .txt file from a squarewave excitation measurement to get an estimate for the QPD sensitivity (beta)
- Lots of changes from Fitting_adam_method_r10.m which you can check out in the code itself
- Can be a bit annoying bc it's set up to work for files which do or don't contain the yPS, but you can set it to always run one way if you always have the same type of file
- RELIES ON LOTS OF OTHER FUNCTIONS NOT FOUND IN THIS FOLDER!
- **requires _oras_codes/Optical Trap Calibration/find_beta_from_sqwv.m_, _loic_codes/adams-method/frf_r9.m_, _loic_codes/adams-method/ps_r6.m_, _loic_codes/adams-method/theor_trap_frf_r9.m_**

##### _sqwv_beta_analysis.m_
- Modified version of Adam's beta analysis code which I can't find anywhere although he claims it exists 
- Takes .txts from squarewave excitation measurement to get an estimate for the QPD sensitivity (beta)
- Same as _oras_codes/Optical Trap Calibration/find_beta_from_sqwv.m_ but just not a function so it's easier to play around with
- Actually doesn't depend on any other functions! How crazy
