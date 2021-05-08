# Toy MC

This project describes a toy MC for estimating the pulse-shape discrimination ability for a wavelength shifting coated component in DEAP-3600.

Requirements: ROOT v5.X

To build:

cd src/

make



--------------------------------
Slow WLS ToyMC Reference Manual:
Author: David Gallacher (dgallacher@snolab.ca)

---------------------------------
To Run:
1.	First run “make” in src/ directory
2.	Then go into scripts/ and run “root make_spectra.C” to ensure PDF histogram files are present in dat/ (make_spectra.C requires rat-deap installation, if input spectra are required email dgallacher@snolab.ca)
3.	To run use: .src/toyMC “seednumber” "type (0==NA,1==NR,2==ER)" “number of trials” “number of experiments” “fileoutname” (Optional) "afterpulsing on/off" (optional)
      * Where you can set the seed for this toy
      * For each “toy” you run there are N experiments of N trials each, specified by the command line inputs
      * Can pass a filename output, defaults to “output.root” and will overwrite in the output/ folder
4.	Run scripts/analyze.C to analyze the output of the toy MC
      * Eg for use: “root -l analyze.C’(“output.root”)’
      * Analyze will look for output files name passed in command line (default is output.root) inside output/, need one for each type
      * Plots from analyze.C are stored in plots/


How it works:

Using gRandom with TRandom3 to handle the random sampling with TH1->GetRandom() and TF1->GetRandom()

* Each toy consists of N experiments, each experiment has N trials
* Inside src/ there are 3 source files; main.cpp, ToyMC.cpp and InverseCDF.cpp, main exists to build ./toyMC , InverseCDF.cpp is a helper class to build inverse CDFs for sampling from, and ToyMC.cpp is where the real work happens
* For each trial we simulate the following:
  * Monoenergetic alphas at 5.3 MeV (Assuming full energy deposition in LAr in the inlet)
  * Apply quenching factor (Set in main.cpp)
  * Convert the energy into “number of photons” by multiplying the edep by the light yield (default = 2k photons/Mev, can be changed) and sampling from a gaussian with mean = edep*Lightyield , and sigma = sqrt(mean)
  * Sample from a uniform distribution for “shadowFraction”, this determines how much of the light gets shadowed, this can vary between 50 and 100% of light shadowed
  * Number of LAr photons = sampledPhotons*(1-shadowFraction)
    * Each photon gets a LAr wavelength sampled from the input spectrum
    * Each photon gets a delta-t sampled from LAr input pulse-shape
    * Add a delta-t from the TPB Pulseshape
    * Change Wavelength to a randomly sampled TPB wavelength
    * If the photon wavelength is < 400 nm, kill it to simulate acrylic absorption
    * Evaluate probability from PMT efficiency curve, throw a uniform number, if number> pmt efficiency then kill photon
    * Pushback final time to “times” vector
    * Evaluate from AP probability distribution, if we have an afterpulse, then pushback the initial pmt hit time + ap time into the “times” array.

  * Number of Pyrene photons = sampledPhotons*shadowFraction
    * Sample from a uniform distribution, if the Random number is > pyrene wavelength shifting efficiency then kill it.
    * Sample from a binomial with mean =2,prob = 0.5; for the number of “reflections”
    * Raise the reflection probability by the number of reflections, sample from a uniform distribution and if the random number > total reflection prob then kill photon
    * Add a delta-t from the pyrenePS Pulseshape
    * Change wavelength to random wavelength from pyrene wavelength distribution
    * Evaluate the acrylic attenuation PDF, this assumes 100cm of LG transmission for the photon, sample a uniform number, if this number > survival prob from acrylic attenuation then kill photon
    * Evaluate probability from PMT efficiency curve, throw a uniform number, if number> pmt efficiency then kill photon
    * If photon survives PMT efficiency then add delta-t from PMT timing distribution
    * Pushback final time to “times” vector
    * Evaluate from AP probability distribution, if we have an afterpulse, the pushback the initial pmt hit time+ ap time into the “times” array.

* Fill TTree with results for each trial.
* TTree Contains:
    * Vector of times of hit photons
    * Number of total hits
    * Number of AP hits
    * Shadow fraction, used to vary “solid angle” effect
    * Number of photons created (Pyrene +LAr)
    * Number of Pyrene photons
    * Trial Number
    * Experiment Number
    * Energy deposited
    * Seed number

We also simulate "Bulk LAr" events in a similar manner for comparison to the inlet-alphas with PPS coating

For Nuclear recoils (type==1) the simulations steps are the same as above, with the following changes:

* Energy deposition sampled from a uniform distribution between 50 and 250 keV
* Apply Ar40 NR Quenching factor
* No shadowing calculated so no "pyrene" photons are produced
* Use pure LAr pulse-shape with NR parameters from rat-deap

For Electronic recoils (type==2) the simulations steps are the same as above, with the following changes:

* Energy deposition sampled from the beta-decay spectrum of Ar39 in dat/spectras.root
* Apply ER Quenching factor (~1.0)
* No shadowing calculated so no "pyrene" photons are produced
* Use pure LAr pulse-shape with ER parameters as implemented in rat-deap and shown in https://arxiv.org/abs/2001.09855
