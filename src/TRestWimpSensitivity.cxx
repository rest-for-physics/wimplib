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

//////////////////////////////////////////////////////////////////////////
/// TRestWimpSensitivity is the main class used to generate senstitivity
/// prospects for WIMP searches experiments. It requires different inputs
/// from an rml file to define different paramented such as the target
/// elements, WIMP density, velocity distributions of the WIMPs in the
/// lab reference frame, exposure, background level, energy spectra for
/// the nuclear recoil, energy range for the sensitivity derivation and
/// the use of the quenching factor.
///
/// #RML definition
///
/// - *addElement* : Add an element to the target material, at least one
/// element should be present with the following parameters:
/// - *nucleusName* : Name of the nucleus to add
/// - *anum* : Atomic number of the nucleus
/// - *znum* : Number of protons in the nucleus
/// - *abundance* : Mass percentage of the nucleus in the target material.
/// - *wimpDensity* : WIMP density in the DM halo in GeV/cm3
/// - *labVelocity* : WIMP velocity in the lab (Earth) frame reference
/// in km/s
/// - *escapeVelocity* : WIMP escape velocity (km/s)
/// - *rmsVelocity* : WIMP RMS velocity (km/s)
/// - *exposure* : Detector exposure in kg*day
/// - *backgroundLevel* : Detector background level in c/keV/day
/// - *energySpectra* : Energy range for the recoil spectra in keV
/// - *energySpectraStep* : Energy step for the recoil spectra in keV
/// - *energyRange* : Energy range for the sensitivity prospects in keV
/// - *useQuenchingFactor* : Use or not quenching factor
/// The following example illustrates the definition of the common WIMP
/// parameters for a Ne+Isobutane mixture at 98% in volume.
///
/// \code
///    <TRestWimpSensitivity name="WIMPSensitivity" title="WIMP Sensitivity" verboseLevel="info">
///        <addElement nucleusName="Ne" anum="20.1797" znum="10" abundance="0.95" />
///        <addElement nucleusName="C" anum="12.0107" znum="6" abundance="0.04" />
///        <addElement nucleusName="H" anum="1.00784" znum="1" abundance="0.01" />
///        <parameter name="wimpDensity" value="0.3" />
///        <parameter name="labVelocity" value="232"/>
///        <parameter name="rmsVelocity" value="220"/>
///        <parameter name="escapeVelocity" value="544"/>
///        <parameter name="exposure" value="116.8"/>
///        <parameter name="background" value="1"/>
///        <parameter name="energySpectra" value="(0,2)"/>
///        <parameter name="energySpectraStep" value="0.01"/>
///        <parameter name="energyRange" value="(0.1,1.1)"/>
///        <parameter name="useQuenchingFactor" value="true"/>
///    </TRestWimpSensitivity>
/// \endcode
///
/// Besides target material elements, the other parameters are optional,
/// and if not provided they will take their default values.
///
/// /// ### Using this class
///
/// Once we have created an instance of this class we will be able to access the
/// different parameters to derive the WIMP sensitivity.
///
/// For example, the following code will perform the derivation of the WIMP
/// sensitivity for a given rml file.
///
/// \code
/// TRestWimpSensitivity WS("wimp.rml");
/// WS.PrintMetadata();
/// double wimpMass = 0.25; //GeV
/// const double sensitivity = WS.GetSensitivity(wimpMass);
///   if(sensitivity > 0 ){
///     std::cout<<"WIMP mass "<<wimpMass<<" Sensitivity "<<sens<< " cm-2"<<std::endl;
///   } else {
///     std::cout<<"Cannot provide sensititivy with this configuration try to
///     incresase the WIMP mass"<<std::endl;
///   }
/// \endcode
///
/// This will return the expected sensitivity for the rml file used above for
/// a WIMP mass of 0.25 GeV.
///
///--------------------------------------------------------------------------
///
/// REST-for-Physics - Software for Rare Event Searches Toolkit
///
/// History of developments:
///
/// 2022-October First implementation
/// JuanAn Garcia
///
/// \class TRestWimpSensitivity
/// \author: JuanAn Garcia   e-mail: juanangp@unizar.es
///
/// <hr>
///

#include "TRestWimpSensitivity.h"
#include "TRestWimpUtils.h"

ClassImp(TRestWimpSensitivity);

using namespace std;

///////////////////////////////////////////////
/// \brief Constructor loading data from a config file
///
/// \param configFilename A const char* giving the path to an RML file.
/// \param name The name of the specific metadata. It will be used to find the
/// corresponding TRestWimpSensitivity section inside the RML.
///
TRestWimpSensitivity::TRestWimpSensitivity(const char* configFilename, const string& name)
    : TRestMetadata(configFilename) {
    Initialize();

    LoadConfigFromFile(fConfigFileName, name);
}
///////////////////////////////////////////////
/// \brief Default destructor
///
TRestWimpSensitivity::~TRestWimpSensitivity() {}

///////////////////////////////////////////////
/// \brief Initialization of TRestWimpSensitivity members
///
void TRestWimpSensitivity::Initialize() {
    SetSectionName(this->ClassName());
    SetLibraryVersion(LIBRARY_VERSION);
}

///////////////////////////////////////////////
/// \brief Initialization of TRestWimpSensitivity
/// members through a RML file
///
void TRestWimpSensitivity::InitFromConfigFile() {
    this->Initialize();
    TRestMetadata::InitFromConfigFile();
    ReadNuclei();
}

///////////////////////////////////////////////
/// \brief Initialization of TRestWimpSensitivity
/// members through a RML file
///
void TRestWimpSensitivity::ReadNuclei() {
    fNuclei.clear();

    TiXmlElement* ElementDef = GetElement("addElement");
    while (ElementDef) {
        TRestWimpNucleus nucleus;
        nucleus.fNucleusName = GetFieldValue("nucleusName", ElementDef);
        nucleus.fAnum = StringToDouble(GetFieldValue("anum", ElementDef));
        nucleus.fZnum = StringToInteger(GetFieldValue("znum", ElementDef));
        nucleus.fAbundance = StringToDouble(GetFieldValue("abundance", ElementDef));
        fNuclei.emplace_back(nucleus);
        ElementDef = GetNextElement(ElementDef);
    }
}

///////////////////////////////////////////////
/// \brief Get recoil spectra for a given WIMP
/// mass and cross section
/// Better performance version
///
std::map<std::string, TH1D*> TRestWimpSensitivity::GetRecoilSpectra(const double wimpMass,
                                                                    const double crossSection) {
    std::map<std::string, TH1D*> recoilMap;

    const int nBins = (fEnergySpectra.Y() - fEnergySpectra.X()) / fEnergySpectraStep;

    for (auto& nucl : fNuclei) {
        std::string histName = "RecoilSpc_" + std::string(nucl.fNucleusName);
        TH1D* recoilSpc =
            new TH1D(histName.c_str(), histName.c_str(), nBins, fEnergySpectra.X(), fEnergySpectra.Y());

        // Build vector of tuples=(recoilEnergy, minimum velocity, rate) used in further calculations
        std::vector<std::tuple<double, double, double>> tEnergyVminRate;
        for (int i = 0; i < recoilSpc->GetNbinsX(); i++) {
            double E = recoilSpc->GetBinCenter(i);
            if (E <= 0) continue;
            tEnergyVminRate.push_back(
                std::make_tuple(E, TRestWimpUtils::GetVMin(wimpMass, nucl.fAnum, E), 0));
        }

        const double nNuclei =
            nucl.fAbundance * TRestWimpUtils::N_AVOGADRO * 1E3 / nucl.fAnum;  // Number of atoms
        const double vMin = std::get<1>(
            tEnergyVminRate.at(0));  // element 0 should be the lowest (positive) energy -> lowest vMin
        const double vMax = fEscapeVelocity + fLabVelocity;

        // calculation of the rate for each recoil energy
        double rate{0};  // will contain integral from minimun vMin to vMax, idem  integral_min(vMin)^vMax
        const double velStep = 0.1;  // km/s
        int j = 0;
        for (double v = vMin; v < vMax; v += velStep) {
            const double flux = 1E5 * v * fWimpDensity / wimpMass;
            const double diffRate =
                flux *
                TRestWimpUtils::GetDifferentialCrossSectionNoHelmFormFactor(wimpMass, crossSection, v,
                                                                            nucl.fAnum) *
                TRestWimpUtils::GetVelocityDistribution(v, fLabVelocity, fRmsVelocity, fEscapeVelocity);
            rate += diffRate * velStep;
            // save (in 3rd element of tEnergyVminRate tuples) the integral from minimun vMin to each vMin,
            // idem integral_min(vMin)^vMin
            while (j < (int)tEnergyVminRate.size()) {
                const double vmin = std::get<1>(tEnergyVminRate.at(j));
                if (vmin < v) {
                    // std::get<2>(tEnergyVminRate.at(j)) = rate; //les precise
                    std::get<2>(tEnergyVminRate.at(j)) = rate - diffRate * (v - vmin);  // more precise
                    j++;
                } else
                    break;
            }
        }
        /*obtain the rate (integral from each vMin to vMax) by substracting integral from minimun vMin to each
           vMin to the integral from minimun vMin to vMax
                idem: integral_vMin^vMax = integral_min(vMin)^vMax - integral_min(vMin)^vMin */
        for (auto& [E, vmin, r] : tEnergyVminRate) {
            if (vmin > vMax) continue;  // r=0
            const double formFactor = TRestWimpUtils::GetHelmFormFactor(E, nucl.fAnum);
            r = (rate - r) * formFactor * formFactor * TRestWimpUtils::SECONDS_PER_DAY * nNuclei;
        }

        // copy results to recoilMap
        j = 0;
        for (int i = 0; i < recoilSpc->GetNbinsX(); i++) {
            const double recoilEnergy = recoilSpc->GetBinCenter(i);
            // const double recoilRate = std::get<2> (tEnergyVminRate.at(i));
            while (j < (int)tEnergyVminRate.size()) {
                if (recoilEnergy == std::get<0>(tEnergyVminRate.at(j))) {
                    recoilSpc->SetBinContent(i, std::get<2>(tEnergyVminRate.at(j)));
                    j++;
                } else
                    break;
            }
        }

        recoilMap[std::string(nucl.fNucleusName)] = recoilSpc;
    }
    return recoilMap;
}

///////////////////////////////////////////////
/// \brief Get sensitivity for a give WIMP mass
///
const Double_t TRestWimpSensitivity::GetSensitivity(const double wimpMass) {
    if (fNuclei.empty()) {
        RESTError << "Please add at least one element to the rml file" << RESTendl;
    }

    if (fUseQuenchingFactor) CalculateQuenchingFactor();

    //const int nBins = (fEnergySpectra.Y() - fEnergySpectra.X()) / fEnergySpectraStep;
    double nMeas = 0;

    const double crossSection = 1E-45;
    auto rSpc = GetRecoilSpectra(wimpMass, crossSection);

    for (auto& nucl : fNuclei) {
        auto recoilSpc = rSpc[std::string(nucl.fNucleusName)];

        for (int i = 1; i < recoilSpc->GetNbinsX(); i++) {
            double recoilEnergy = recoilSpc->GetBinCenter(i);
            const double recoilRate = recoilSpc->GetBinContent(i);

            if (fUseQuenchingFactor)
                recoilEnergy *= quenchingFactor[std::string(nucl.fNucleusName)]->GetBinContent(i);

            if (recoilEnergy < fEnergyRange.X() || recoilEnergy > fEnergyRange.Y()) continue;
            nMeas += recoilRate * fEnergySpectraStep;
        }
    }

    double bckCounts = 0;

    auto recSpc = rSpc[std::string(fNuclei.front().fNucleusName)];
    for (int i = 1; i < recSpc->GetNbinsX(); i++) {
        const double en = recSpc->GetBinCenter(i);
        if (en < fEnergyRange.X() || en > fEnergyRange.Y()) continue;
        bckCounts += fBackground * fEnergySpectraStep;
    }
    bckCounts *= fExposure;

    for (auto& [name, histo] : rSpc) delete histo;
    rSpc.clear();

    if (nMeas == 0) return 0;

    double signalCounts = 0, prob = 0;

    do {
        prob = 0;
        for (int mu = signalCounts; mu < (signalCounts + bckCounts + 10000); mu++) {
            if (bckCounts <= 1.e3)
                prob += TMath::Poisson(mu + bckCounts, bckCounts);
            else if (bckCounts > 1.e3)
                prob += TMath::Gaus(mu + bckCounts, bckCounts, TMath::Sqrt(bckCounts), true);
        }
        signalCounts++;
    } while (fabs(prob - 0.1) > 0.01 && signalCounts < 1E6);

    const double sensitivity = signalCounts * 1E-45 / (nMeas * fExposure);

    return sensitivity;
}

///////////////////////////////////////////////
/// \brief Calculate Quenching factor and
/// stores in a map
///
void TRestWimpSensitivity::CalculateQuenchingFactor() {
    if (!quenchingFactor.empty()) return;

    std::cout << "Calculating quenching factor " << std::endl;

    const int nBins = (fEnergySpectra.Y() - fEnergySpectra.X()) / fEnergySpectraStep;

    for (auto& nucl : fNuclei) {
        std::string histName = "QF_" + std::string(nucl.fNucleusName);
        TH1D* QF =
            new TH1D(histName.c_str(), histName.c_str(), nBins, fEnergySpectra.X(), fEnergySpectra.Y());
        for (int i = 1; i < QF->GetNbinsX(); i++) {
            const double recoilEnergy = QF->GetBinCenter(i);
            const double qF = TRestWimpUtils::GetQuenchingFactor(recoilEnergy, nucl.fAnum, nucl.fZnum);
            if (!isnan(qF) && qF > 0)
                QF->SetBinContent(i, 1. / qF);
            else
                QF->SetBinContent(i, 0);
        }
        quenchingFactor[std::string(nucl.fNucleusName)] = QF;
    }
}

///////////////////////////////////////////////
/// \brief Return output file format with different
/// parameters used in the calculation.
///
const std::string TRestWimpSensitivity::BuildOutputFileName(const std::string& extension) {
    std::stringstream ss;
    ss << "WimpSensitivity_";

    for (auto& nucl : fNuclei) ss << nucl.fNucleusName << "_" << nucl.fAbundance << "_";

    ss << "WD_" << fWimpDensity << "_";
    ss << "Vel_" << fLabVelocity << "_" << fRmsVelocity << "_" << fEscapeVelocity << "_";
    ss << "Bck_" << fBackground << "_";
    ss << "Exp_" << fExposure << "_";
    ss << "RecEn_" << fEnergySpectra.X() << "_" << fEnergySpectra.Y() << "_" << fEnergySpectraStep << "_";
    ss << "EnRange_" << fEnergyRange.X() << "_" << fEnergyRange.Y() << "_";

    if (fUseQuenchingFactor)
        ss << "usingQF";
    else
        ss << "noQF";

    ss << extension;

    std::cout << "Output File " << ss.str() << std::endl;

    return ss.str();
}

///////////////////////////////////////////////
/// \brief Return a map of histograms with the
/// Form Factor of the different elements.
///
std::map<std::string, TH1D*> TRestWimpSensitivity::GetFormFactor() {
    std::map<std::string, TH1D*> formFactor;

    const int nBins = (fEnergySpectra.Y() - fEnergySpectra.X()) / fEnergySpectraStep;

    for (auto& nucl : fNuclei) {
        std::string histName = "FF_" + std::string(nucl.fNucleusName);
        TH1D* FF =
            new TH1D(histName.c_str(), histName.c_str(), nBins, fEnergySpectra.X(), fEnergySpectra.Y());
        for (int i = 1; i < FF->GetNbinsX(); i++) {
            const double recoilEnergy = FF->GetBinCenter(i);
            const double helmFF = TRestWimpUtils::GetHelmFormFactor(recoilEnergy, nucl.fAnum);
            FF->SetBinContent(i, helmFF * helmFF);
        }
        formFactor[std::string(nucl.fNucleusName)] = FF;
    }

    return formFactor;
}

///////////////////////////////////////////////
/// \brief Prints on screen the details about WIMP
/// parameters, stored in TRestWimpSensitivity.
///
void TRestWimpSensitivity::PrintMetadata() {
    TRestMetadata::PrintMetadata();

    for (auto& nucl : fNuclei) nucl.PrintNucleus();

    RESTMetadata << "WimpDensity: " << fWimpDensity << " GeV/cm3" << RESTendl;
    RESTMetadata << "WimpVelocity; VLab: " << fLabVelocity << " VRMS: " << fRmsVelocity
                 << " VEscape: " << fEscapeVelocity << " km/s" << RESTendl;
    RESTMetadata << "Exposure: " << fExposure << " kg*day" << RESTendl;
    RESTMetadata << "Background Level: " << fBackground << " c/keV/day" << RESTendl;
    RESTMetadata << "Recoil energy range: (" << fEnergySpectra.X() << ", " << fEnergySpectra.Y()
                 << ") Step: " << fEnergySpectraStep << " keV" << RESTendl;
    RESTMetadata << "Sensitivity energy range: (" << fEnergyRange.X() << ", " << fEnergyRange.Y() << ") keV"
                 << RESTendl;
    RESTMetadata << "Use quenching factor: " << fUseQuenchingFactor << RESTendl;
    RESTMetadata << "+++++" << RESTendl;
}
