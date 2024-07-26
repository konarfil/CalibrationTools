#ifndef CALIB_INFO_H
#define CALIB_INFO_H

struct calib_info
{
  double a; // first calibration parameter
  double b; // second calibration parameter
  double chi2_NDF_A; // chi-squared over NDF of the fit of the 482 keV peak for the best found a and b
  double chi2_NDF_B; // chi-squared over NDF of the fit of the 976 keV peak for the best found a and b
  double loss; // value of the loss function for the best found a and be
};
#endif