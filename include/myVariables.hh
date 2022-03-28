#include <fstream>
#include <vector>
#include <iostream>
#include <string>

void readConfigFile();
void readAlphaFile();
void readGammaFile();
void readXRayFiles();

std::vector<double> alphaEnergy, alphaIntensity, gammaEnergy, kEnergy, kIntensity;
std::vector< std::vector<double> > gammaIntensity, ceEnergy, ceIntensity, lEnergy, lIntensity;

double annularPositionZ, circularPositionZ;
int Z, optionAlpha, optionGamma, optionCE, optionX, nAlphas, nGammas;