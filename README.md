# emily-codes
Emily's codes for live cell organelle motility, OT, and immunofluorescence analysis.

Note: most of these codes reference the functions in the General_codes folder. Ensure you add all of these functions to a single folder and addpath to this directory before trying to run your analyses.

Organelle_motility_TrackMate

	1. Analyze_Tracjmate_motility_v2_em_master
	Adapted from Abdullah's code. Takes in a folder of .csv files from the Spots_in_tracks_statistics files in TrackMate's results and generates a .mat file with the relevant information for future analysis. The results must have the 'cell centroid' defined or the x and y coordinates of the axon location closest to the cell soma in the rightmost columns. Additionally, the axon length generated when making the kymograph should be added to the column after the centroid, so that if needed the data can be normalized to axon length. 

	2. Optimal_multi_min_lr_per_cell
	COMPUTATIONALLY EXPENSIVE, takes long to run (use analysis computer). Takes in the .mat files from step 1 and iterates through multiple potential minimum run lengths for analysis. It plots the 	diffusive and processive alphas for each minimum run length as well as the directional bias at each set minimum run length within a range. This helps determine the appropriate minimum run length for a 	processive run to be used in step 3. Note, you can still run step 3 and 4 with a guess on the minimum 	run length before running this code, but this code will make the guess more accurate and justified. Used in both Emily's manuscripts in supplementary figures. 

	3. Analyze_dispersion_Rg_MSD_per_cell_5cond_v2
	Adapted from Abdullah's code. Takes in the .mat files from step 1 and generates .mat files of the Rg, MSD, and directional bias per cell to be read by the next code. The _kymo_coords version of this code uses a version of the data that has been flattened along the axis of the kymograph and performs the same analysis. 

	4.Plot_Rg_and_MSD_analysis_per_cell_4cond_w_ax_len_bdnf_v5
	Adapted from Abdullah's code. Takes in the .mat files from step 3 and generates plots of the alpha, radius of gyration, and directional bias. It can also normalize the data to axon length if that was included in step 1. It plots the diffusive alpha to show that the set minimum run length analysis is appropriate.

	Plot_xk_yk_per_cell_all_em_updated
	Adapted from Abdullah's code. This plots all the trajectories from the analysis on one plot in the positions from the original data in the .mat file from step 1. 

	Plot_xk_yk_per_cell_stationary_track
	Adapted from Abdullah's code. This isolates the stationary trajectories by their displacements and their alpha values and plots them. The distribution fit is then used to calculate the tracking uncertainty for a 	cargo based on the stationary data from TrackMate. Used in Emily's polyQ huntingtin manuscript. 

	Plot_trajectories_color_gradient_v2_em_all
	Adapted from Abdullah's code. This plots all the trajectories from the analysis on one plot, normalized around zero using the .mat files from step 1. Used in Emily's S421D paper Figure 2. 

	get_axon_lengths_per_cell
	Generates a .mat file with the axon lengths from each of the .csv files, to be used in the flux analysis. 

	TrackMate_to_kymo_and_flux_v4
	Takes in the .mat files from step 1 and analyzes the flux 5 boxes of a specified box_width (variable is half the width of the box you want). Plots the resulting fluxes and saves a variable with these fluxes that can be used to make nicer plots (flux_inset_plot). Note that exit fluxes were the ones used in analysis and this flux counts trajectories that either start inside the box and exit outside the box or those that start outside the box and finish inside it. It also performs a one-way and multi-way ANOVA analysis on the fluxes.

	Analyze_processive_motility_em_bstrp_4kchoose
	Adapted from Abdullah's code. Takes in the .mat files from either TrackMate or KymoButler and plots data per trajectory above a minimum run length. Used mainly to analyze velocity and run length, as well as plots a histogram of their product. Used in Figures 2-6 in Emily's polyQ huntingtin manuscript. 

Organelle_motility_kymobutler



	1. Kymobutler_reader
	Takes in a list of kymobutler files and generates a .mat file with the relevant parameters to be used in the next codes. Note: you will need to have run the get_axon_lengths_file to use the flux analysis (see TrackMate codes). 

	2. TrackMate_to_kymo_and_flux_v4_KB
	Takes in the .mat files from kymobutler_reader and analyzes the flux 5 boxes of a specified box_width (variable is half the width of the box you want). Plots the resulting fluxes and saves a variable with these fluxes that can be used to make nicer plots (flux_inset_plot). Note that exit fluxes were the ones used in analysis and this flux counts trajectories that either start inside the box and exit outside the box or those that start outside the box and finish inside it. It also performs a one-way and multi-way ANOVA analysis on the fluxes. This code was used in Emily's polyQ huntingtin manuscript. You can also modify this code to run trackmate files, you just need to uncomment the section that turns the 2D positions to 1D. 
	2-alternative. TrackMate_to_kymo_and_flux_v5_KB_alt_more_stringent
	This code is as the previous one but has a more stringent flux definition that only counts trajectories that move from one side to the other side of the box as a directional flux.
	3. flux_inset_plot
	This code makes the plots from step 2 using the .mat files generated from step 2 so that you don't need to iterate through all the files again. 
	4. flux_box_size_trends
	This code is used to show the trends in flux when the box width is modified. Running step 2 for multiple box widths and including the box width in the .mat filename is required before running this code. This code then takes the results from step 2 and plots the trends in relevant parameters as box widths increase. Used in Figures S2-6 in Emily's polyQ huntingtin manuscript. 
	5. TrackMate_to_kymo_and_colour_coded_flux_KB
	This code is used to generate the plots of flux analysis itself with the colours representing flux direction and the boxes on top. It is a great resource for testing your flux analysis parameters and spot checking results. It runs on only the sample figure files because otherwise your computer will be angry from all of the iterations and plots being generated. Used for Figures 2-6 in Emily's polyQ huntingtin manuscript.

	Analyze_processive_motility_em_bstrp_4kchoose
	Adapted from Abdullah's code. Takes in the .mat files from either TrackMate or KymoButler and plots data per trajectory above a minimum run length. Used mainly to analyze velocity and run length, as well as plots a histogram of their product. Used in Figures 2-6 in Emily's polyQ huntingtin manuscript. 

	Analyze_dispersion_Rg_MSD_per_cell_KB_v1
	Adapted from Abdullah's code. Takes in the .mat files from the reader and generates files for plotting the alpha and radius of gyration in the Plot_Rg_and_MSD... code. Not used in published work. 

OT
	Calibration
	
	Most of these codes are Ora's codes just modified in Rbead for Nanodiamonds. See Ora's ReadMe for detailed instructions. 

		Plotting_raw_AOD
		Adapted from Magda's codes. Plots the raw AOD signal in both on and off axis directions.

		Optimizing_PS_fits
		This code takes in a list of good calibrations and attempts to find the best match to a calibration you're having trouble fitting because of the magnitude response. It measures the residuals from the raw power spectrum (blue) compared to the power spectra from the traces you got good fits from (plots this in magenta and the region it is using for the residual calculation in yellow) and gives the index of the file that fits best to be used as the initial guess in the fitting code.

	Force_traces_analysis

	1. force_trace_analysis_fiona_v5_new_em_neuron
	Adapted from Abdullah's code. Runs one trace at a time. Takes in a force trace and the beta and ktrap from its calibration (entered manually), as well as the direction of the microtubule it is on (x, y, or xy). It then calculates the stall force distribution and plots the force trace with forces and or distances and generates a .mat file that can be used in future analysis. See the 	Steps_to_calibrate_and_run_OT_txt_files for information on when to flip plot or not. 

	1-alt. force_trace_analysis_fiona_r6_em_neuron
	Adapted from Adam's code. Does the same thing as step 1 but makes a slightly different variable in the end and doesn't use the flip plot concept. Instead the r_cam_ant variable must be used to define microtubule orientation, and you must manually select the zero in both x and y, then select the microtubule direction from perinuclear to peripheral on a scatter of all x points. 

	2. force_trace_cumulative_r6_em_3_neuron
	Adapted from Abdullah's code. Plots the force distribution from all of step 1's .mat files and attempts to fit multiple gaussians to it. Is set to compare 2 conditions. 

	3. estimate_unbinding_rate_em_a_3
	Adapted from Adam's codes. Uses the Berger lab paper from 2019 to estimate unbinding rates from the force trace .mat files from step 1. Emily's version references the fit_single_exponential to fit a single exponential to each half of the data (above and below zero). 

	4. estimate_binding_and_unbinding_rate_motors_em_3_neuron
	Adapted from Abdullah's version of Adam's codes. Plots binding and unbinding rates from the .mat files in step 1, but I usually use the binding rate figures from this code and the unbinding rate data from step 3. 

	The OT force trace analysis codes are currently not very efficient, it would be nice in the future to consolidate the codes from steps 2-4 but I haven't had the time to do that.


IF_analysis_plotting 

	Cell_morphology_plots
	Takes a list of .csv files generated by ImageJ with various parameters from the analyze particles analysis and plots histograms and violin plots of the distribution from each cell. Used in Figure S1 in Emily's S421D huntingtin paper.

	intensity_profile_three_colour
	Takes in a list of .csv files generated from an intensity linescan of a 3-channel image in imageJ and plots the linescans.

	motor_colocalization_2
	Takes in a folder of single channel images and another folder of single channel images made from multi-colour experiments. Analyzes intensity colocalization for each channel. Used in Figure 4 and Figure S5 of Emily's S421D huntingtin paper. Can be modified to measure emission crosstalk.

	Htt_tubb3_int_arb_analysis
	Takes in .csv files with intensities or arborization lengths in neurons and plots them. Also performs a one-way ANOVA on the results. Used in Figure 1 in Emily's polyQ huntingtin manuscript.
