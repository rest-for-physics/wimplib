/*************************************************************************
 * This file is part of the REST software framework.                     *
 *                                                                       *
 * Copyright (C) 2016 GIFNA/TREX (University of Zaragoza)                *
 * For more information see http://gifna.unizar.es/trex                  *
 *                                                                       *
 * REST is free software: you can redistribute it and/or modify          *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation, either version 3 of the License, or     *
 * (at your option) any later version.                                   *
 *                                                                       *
 * REST is distributed in the hope that it will be useful,               *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          *
 * GNU General Public License for more details.                          *
 *                                                                       *
 * You should have a copy of the GNU General Public License along with   *
 * REST in $REST_PATH/LICENSE.                                           *
 * If not, see http://www.gnu.org/licenses/.                             *
 * For the list of contributors see $REST_PATH/CREDITS.                  *
 *************************************************************************/

#ifndef RestCore_TRestWimpSensitivity
#define RestCore_TRestWimpSensitivity

#include <TH1D.h>
#include <TRestMetadata.h>
#include <TRestWimpNucleus.h>

/// Container class for WIMP metadata
///
class TRestWimpSensitivity : public TRestMetadata {
   private:
    /// A vector containing TRestWimpNucleus with different nucleus parameters
    std::vector<TRestWimpNucleus> fNuclei;
    /// WIMP density in GeV/cm3
    Double_t fWimpDensity = 0.3;
    /// WIMP velocity in the lab (Earth) frame reference in km/s
    Double_t fLabVelocity = 232;
    /// WIMP escape velocity (km/s)
    Double_t fEscapeVelocity = 544;
    /// WIMP RMS velocity (km/s)
    Double_t fRmsVelocity = 220;
    /// Detector exposure in kg*day
    Double_t fExposure = 365. * 0.32;
    /// Detector background level in c/keV/day
    Double_t fBackground = 1;
    // TODO add option to use a histogram for the spectra
    /// Energy range for the recoil spectra in keV
    TVector2 fEnergySpectra = TVector2(0, 2);
    /// Energy step for the recoil spectra in keV
    Double_t fEnergySpectraStep = 0.01;
    /// Energy range for the sensitivity prospects in keV
    TVector2 fEnergyRange = TVector2(0.1, 1.1);
    /// Use or not quenching factor
    Bool_t fUseQuenchingFactor = true;

    /// Map containing quenching factor for a nucleus
    std::map<std::string, TH1D*> quenchingFactor;  //!

   public:
    TRestWimpSensitivity(const char* configFilename, const std::string& name = "");

    ~TRestWimpSensitivity();

    void Initialize() override;
    void InitFromConfigFile() override;
    void PrintMetadata() override;

    void ReadNuclei();
    const Double_t GetSensitivity(const double wimpMass);
    void CalculateQuenchingFactor();
    const std::string BuildOutputFileName(const std::string& extension = ".txt");

    std::map<std::string, TH1D*> GetRecoilSpectra(const double wimpMass, const double crossSection);
    std::map<std::string, TH1D*> GetFormFactor();
    inline auto GetQuenchingFactor() { return quenchingFactor; }

    ClassDefOverride(TRestWimpSensitivity, 1);
};

#endif
