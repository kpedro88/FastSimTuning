Generation:
	SinglePi_SIM_split_temp.py is a template python configuration file for generating full simulation pion samples.
	FStempEsplit.sh and FSsubsplit.sh are shell scripts which create actual configuration files from the template for specific energies.
	(These files are "split" for parallel batch submission of jobs.)

Analysis:
	fullsimpionanalyzer_split*_temp_cfg.py are template python configuration files for several different splittings.
	FAtempsplit.sh and FArun.sh are shell scripts which create actual configuration files from the template for specific energies.
	These config files run the analyzer FullSimPionAnalyzer, which computes the reconstructed energy for each event in ECAL and HCAL using PCaloHits.
	The HCAL sampling factors come from SimCalorimetry/HcalSimProducers/python/hcalSimParameters_cfi.py.
	(The calculation of the sampling factors is described at https://twiki.cern.ch/twiki/bin/viewauth/CMS/HcalSamplingFactors.)
	See https://twiki.cern.ch/twiki/bin/viewauth/CMS/HcalSignalEvaluation for a motivation of the Poisson smearing process used to reconstruct energy in HF.
	
	fullsimetatimeanalyzer_50_split5000_cfg.py is a python configuration file which runs the analyzer FullSimEtaTimeAnalyzer on the 50 GeV pion sample.
	This analyzer creates histograms of the PCaloHit timing distribution for each ieta in HCAL, with each entry weighted by the PCaloHit energy.

Macros:
	fs_energy_res_cballD_lim.C shows the basic routine for fitting energy distributions and storing the fit parameters.
	fs_energy_res_cballD_lim_custom.C is the same macro with "customizations" added to nudge fits away from bad local minima.
	fs_energy_res_cballD_comp.C is a macro which allows the comparison of fit parameters and results for a given ieta bin, or for all ieta bins.
	fs_energy_res_cballD_interp_safe.C is the same as the previous macro, but with an interpolation routine to generate parameters for intermediate energies.
	
	get_time_peaks.C finds the peaks of the timing histograms created by FullSimEtaTimeAnalyzer and prints them in python format.

All macros should be loaded and compiled before use:
root -l
.L [macro].C+