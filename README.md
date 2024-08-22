# CalibrationTools
 The repository contains tools for energy calibration of the SuperNEMO detector. The calibration procedure is divided into three parts: Falaise modules CalibrationCutsModule and Charge2EnergyModule, and command line application sncalib. CalibrationCutsModule and sncalib are used to calculate calibration parameters based on measured calibration data, while Charge2EnergyModule uses the calibration parameters to transform measured charge into energy and save it into the CD bank. The repository also contains the EnergyCorrectionCalculator class which implements both optical correction and energy loss correction. These corrections are then considered in the energy calibration. For details about energy corrections and the whole calibration process see my master thesis at docDB #5933.

## Instalation (CC-IN2P3)
First load the current Falaise enviroment. At the time of writing this README, it is done by the following commands.
```
source /sps/nemo/sw/snswmgr/snswmgr.conf
snswmgr_load_setup falaise@5.1.2
```

With loaded enviroment the calibration tools are installed by running this set of commands:

```
git clone https://github.com/konarfil/CalibrationScript.git
cd CalibrationScript
mkdir build
cd build
cmake ..
make
```

## CalibrationCutsModule
CalibrationCutsModule reads a brio file with calibration data, applies data cuts, extracts information necessary for the calibration and saves it into a root file. The input brio file must contain the pCD bank with measured charge, the PTD bank with track information and the CD bank for track-calo associations. At the current version the module uses two main cut conditions to choose tracks for the calibration: the source vertex must be close to a calibration source and the track must have exactly one associated calorimeter hit. The module also skips kinked tracks. These can be maybe included in the future, but the energy loss correction might not work for them that well. **It is important to ensure that the calibration source positions match with the plasma propagation model used for track reconstruction. It is possible to either load calibration source positions from Falaise (these curently do not match the plasma propagation model) or load them from a text file by adding the file path to the configuration file ("source_pos_path" variable).**

The output root file contains a TTree for each optical module (OM). Each entry in a TTree corresponds to one track and contains measured charge, source vertex position, calorimeter vertex position and calorimeter vertex position in coordinate system of the OM (this is later used for optical correction). This is minimal information needed for the calibration but the list of variables can be expanded in the future to also serve diagnostic purposes.

### Usage
`flreconstruct -i /path/to/input.brio -p /path/to/config.conf`

See the example configuration file calib_cuts.conf

## SNCalib
SNCalib calculates calibration parameters of individual OMs based on data from the root file produced by CalibrationCutsModule. The calibration parameters are then saved into a csv file. It is expected that the SNCalib program will be used on root files containing larger data (at least 10 hours of measurement). The ROOT hadd command can be used to easily combine data from multiple calibration runs. The calibration algorithm is based on minimization of a loss function of the calibration parameters *a* and *b*. The loss function expresses how far the peaks in the calibration spectrum are from the correct energies for chosen pair of *a* and *b*. We are looking for *a* and *b* for which the loss function is as close to 0 as possible.

At the current version the output csv file contains 6 columns OM_number, a, b, chi2_A, chi2_B, loss. The first three columns contain a pair of calibration parameters for given OM. The other three values are meant for diagnostics. chi2_A and chi2_B are Chi-squared over NDF values of fits of the calibration spectrum. Large chi2_A and chi2_B is typically caused by a defect in the spectrum such as big gain variation. The loss value is the lowest value of the loss function found for given OM. Ideally this value should be lower than the minimization_threshold value set in the config file. A higher value means that the algorithm did not converge in less than max_iterations (also set in the config file). A loss value only slightly higher than minimization_threshold can still mean a decent calibration, although very high values should be investigated.

Besides the three diagnostic values it is also possible to save png images of the fitted calibration spectra by including the "-s" flag (see Usage). These will be saved into the SNCalib/Fits folder. Using the values of chi2_A, chi2_B and loss as well visual examination of the spectra, it is possible to find OMs with unsuccessful calibration.

### Usage
`./sncalib -i /path/to/input.root -o /path/to/output.csv -p /path/to/config.conf -n 90 -s -V`

The program can be used in two modes. If the "-n" argument is included, only the OM with given number (90 in this case) is calibrated and the result is printed. If the "-n" argument is ommited, the program calibrates all modules (that are available in "input.root") and saves the results into "output.csv". Parameters of the tracking gas and the calibration algorithm itself can be set in the "config.conf" file. See the example configuration file "params.conf". The "-s" flag makes the program save pngs of the fitted spectra. The "-V" flag makes the program plot the number of processed OM.

## Charge2EnergyModule
Charge2EnergyModule reads a brio file and based on calibration parameters from a csv file produced by SNCalib transformes measured charges into energies. The brio file must contain pCD bank with measured charge as well PTD and CD banks with information about reconstructed tracks and track-calo associations. Calibration parameters are saved into calibrated_calorimeter_hits inside the CD bank. The energy of the calorimeter hit is set to *aQ + b* where *Q* is the measured charge. This is the energy without any corrections. The corrected energies are saved into auxiliaries of the calorimeter hit as "Ef", "Ef_bc", "Ef_optical" and "Ef_optical_loss" being energies without any corrections, with Birks and Cherenkov correction, with optical correction and with both optical and energy loss correction respectivelly. "Ef_optical_loss" should, thus, provide the best estimate of the initial energy of the particle. If a vertex or a whole associated track is missing than it is impossible to calculate some of the corrections and "Ef_optical" and "Ef_optical_loss" can, thus, also be missing.

### Usage
`flreconstruct -i /path/to/input.brio -p /path/to/config.conf -o /path/to/output.brio`

See the example configuration file charge2energy.conf.

## Contact
In case of any questions or problems, do not hesitate to contact me (Filip Koňařík) through Slack or email (konarfil@cvut.cz).
