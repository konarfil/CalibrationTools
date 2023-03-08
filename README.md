# CalibrationScript
The repository contains the first version of the calibration algorithm. For simplicity we consider the relationship between energy and charge to be $E = kQ$ where $k$ is an unknown calibration parameter and we are trying to find $k$ based on the position of the 976 keV peak(ignoring 482 keV for now). It is divided into 2 scripts - extract_OM_spectra.cxx, and fit_spectra.cpp.
## extract_OM_spectra.cxx
This script selects events containing Bi207 electrons using tracking algorithm and saves energy spectra measured by individual optical modules. Each event is evaluated in 3 steps:
1. Try to find a track on each side of the tracker
2. For each track, check if it meets cut conditions
3. If the track passes all conditions, save energy

Tracking is performed using algorithm developed by Tomáš Křižák (https://github.com/TomasKrizak/TKEvent). It gives us one line for each side of the tracker. From each line we can extract electrons initial vertex (at $x=0$) and its final vertex (at the front face of an OM). Positions of these vertices are than used in the cut.
To pass the cut an electron has to meet several criteria:

- initial vertex within 10 cm from a calibration source
- final vertex within 20 cm from the front face center of the closest triggered OM
- at least 3 triggered cells near the calibration source (checks a rectangle sized 11x10 tracker cell radii)
- less than 16 tracker cells in time coincidence with the OM (time window between -0.2 and 5.0 us)

If the track passes all conditions we calculate its estimated initial energy as
$E_0=kQ+\Delta E(kQ,d;p_1, p_2, p_3, p_4, p_5)=kQ+(p_1d+p_2)e^{p_3kQ}+p_4d+p_5$,

where $\Delta E$ is a correction for energy lost by the electron while passing the tracker. Its form as well as parameters $p_1,...,p_5$ were determined using Falaise simulation. A search for more precise correction model is in progress. Since we do not know the $k$ parameter we calculate the energy for 100 different values of $k$ and as a result save 100 spectra for each OM. From these spectra we than choose the one where the second peak is closest to 976 keV. Its $k$ parameter is the one we are looking for.
## fit_spectra.cpp
A root macro used to fit energy spectra saved by "exctract_OM_spectra" using sum of three gaussian functions. It looks for a fit where the mean value of the first gaussian is closest to 976 keV. This fit is than saved as png image and the $k$ parameter is save into a text file.
## Usage
Following steps describe how to run the scripts on CC-IN2P3:
1. Clone the repository
~~~~
git clone https://github.com/konarfil/CalibrationScript.git
cd CalibrationScript
~~~~
2. Install TKevent
~~~~
$ git clone https://github.com/TomasKrizak/TKEvent.git
$ cd TKEvent/RED_to_TK/
$ source load_environment_new.sh
$ cd ../TKEvent/
$ chmod 755 install.sh
$ ./install.sh
~~~~
3. Compile exctract_OM_spectra
~~~~
$ cd ../../build
$ chmod 755 compile.sh
$ ./compile.sh
~~~~

After this setup you can run the calibration scripts using "run.sh". The bash script asks for a number of run which you want to read and how many events you want to read. Than it creates directories to save the output of extract_OM_spectra and fit_spectra and submits a job to CC-IN2P3. **The script only works for runs 813 and newer !** After the job finishes spectra histograms are saved in "../runs/run_###/OM_histos/", fits are saved in "../runs/run_###/fits/" and values of the $k$ parameter for individual modules are saved in "../runs/run_###/fits/calibration_params.txt". If you log off CC-IN2P3 you should run "load_enviroment.sh"("build" folder) before running "run.sh". 

If you find an error or have any questions feel free to contact me(Filip Koňařík) through Slack or email(konarfil@cvut.cz).
