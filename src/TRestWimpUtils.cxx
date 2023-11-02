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
/// TRestWimpUtils defines several functions to calculate different WIMP
/// parameters, based on WimpRoot code written by I.G. Irastorza, some
/// functions are derived from a D. Diez-Ibañez python code.
///
/// Note units are keV for the recoil energy, GeV for the WIMP mass,
/// km/s for the WIMP velocity and cm-2 for the cross section.
///--------------------------------------------------------------------------
///
/// REST-for-Physics - Software for Rare Event Searches Toolkit
///
/// History of developments:
///
/// 2022-October First implementation
/// JuanAn Garcia
///
/// \class TRestWimpUtils
/// \author: JuanAn Garcia   e-mail: juanangp@unizar.es
///
/// <hr>
///

#include <TMath.h>
#include <TRestWimpUtils.h>

//////////////////////////////////////////////////
/// \brief Get relative nuclear cross section
/// within a WIMP and a nucleon, assuming
/// SCALAR INTERACION
///
const double TRestWimpUtils::GetRelativeNuclearCS(const double wimpMass, const double Anum) {
    const double reducedMass = GetReducedMass(wimpMass, Anum);
    const double reducedMassSingle = GetReducedMass(wimpMass, 1.);
    return (Anum * reducedMass / reducedMassSingle) * (Anum * reducedMass / reducedMassSingle);
}

//////////////////////////////////////////////////
/// \brief Get WIMP-nucleon reduced mass
/// (WIMP mass in GeV)
///
const double TRestWimpUtils::GetReducedMass(const double wimpMass, const double Anum) {
    // WIMP mass in GeV
    return wimpMass * GEV_PER_UMA * Anum / (wimpMass + Anum * GEV_PER_UMA);
}

//////////////////////////////////////////////////
/// \brief Get Helm form factor for a given recoil
/// energy and nucleous target, mass recoil
/// energy in keV
///
const double TRestWimpUtils::GetHelmFormFactor(const double recoilEnergy, const double Anum) {
    // Momentum transfer in keV
    const double q = sqrt(2. * Anum * GEV_PER_UMA * 1E6 * recoilEnergy);
    const double s = 0.9;  // Femtometers-Skin thickness of the nucleus
    const double qs = q * s / HC_KEV_FM;
    // Parametrization of atomic nuclei
    const double RN = sqrt(pow(1.23 * std::cbrt(Anum) - 0.60, 2) +
                           (7. / 3.) * pow(TMath::Pi(), 2) * 0.52 * 0.52 - 5. * s * s);
    const double qR = q * RN / HC_KEV_FM;
    // First Bessel function
    const double bessel1 = sin(qR) / (qR * qR) - cos(qR) / qR;
    // Form factor
    double formFactor = 3. * bessel1 / qR * exp(-0.5 * qs * qs);

    return formFactor;
}

//////////////////////////////////////////////////
/// \brief Get minimum velocity distribution (km/s)
/// for a WIMP to create a recoil energy (keV)
/// higher than recoilEnergy
///
const double TRestWimpUtils::GetVMin(const double wimpMass, const double Anum, const double recoilEnergy) {
    const double reducedMass = GetReducedMass(wimpMass, Anum);
    return sqrt(Anum * GEV_PER_UMA * recoilEnergy * 1E-6 / (2. * reducedMass * reducedMass)) * LIGHT_SPEED;
}

//////////////////////////////////////////////////
/// \brief Get velocity distribution for a given
/// Velocity distribution in Earth frame, f(velocity) in units of 1/(km/s), velocity in km/s
/// Result is the integral for all solid angle of the Boltzmann velocity distribution in Earth frame*/
///
const double TRestWimpUtils::GetVelocityDistribution(const double v, const double vLab, const double vRMS,
                                                     const double vEscape) {
    if (v > vLab + vEscape) return 0;

    const double vAdim = vEscape / vRMS;
    const double Nesc =
        erf(vAdim) - 2. / sqrt(TMath::Pi()) * vAdim *
                         exp(-vAdim * vAdim);  // Nesc=1 for vEscape=infinity (see Lewin&Smith appendix 1a)
    // xMax = max(cosTheta) to meet [vec(v) + vec(vLab)]^2 < vEscape^2 boundary condition (also xMax in
    // [-1,+1])
    const double xMax =
        std::max(-1., std::min(1., (vEscape * vEscape - vLab * vLab - v * v) / (2. * vLab * v)));

    return v / Nesc / (vLab * vRMS * sqrt(TMath::Pi())) *
           (exp(-(v - vLab) * (v - vLab) / (vRMS * vRMS)) -
            exp(-(v * v + vLab * vLab + 2 * v * vLab * xMax) / (vRMS * vRMS)));
}

//////////////////////////////////////////////////
/// \brief Differential cross section without Helm form factor in energy,
/// in [cm/keV]. E in keV, velocity in km/s,
/// wimp mass in Gev/c^2, cross section per
/// nucleon in cm^2, Anum in atomic units (amu)
/// (or atomic weight)
/// Useful function for performance enhancement
///
const double TRestWimpUtils::GetDifferentialCrossSectionNoHelmFormFactor(const double wimpMass,
                                                                         const double crossSection,
                                                                         const double velocity,
                                                                         const double Anum) {
    const double cs = GetRelativeNuclearCS(wimpMass, Anum) * crossSection;
    const double reducedMass = GetReducedMass(wimpMass, Anum);
    const double Emax = 1E6 / LIGHT_SPEED / LIGHT_SPEED * 2. * reducedMass * reducedMass * velocity *
                        velocity / (Anum * GEV_PER_UMA);

    return cs / Emax;
}

//////////////////////////////////////////////////
/// \brief Differential cross section in energy,
/// in [cm/keV]. E in keV, velocity in km/s,
/// wimp mass in Gev/c^2, cross section per
/// nucleon in cm^2, Anum in atomic units (amu)
/// (or atomic weight)
///
const double TRestWimpUtils::GetDifferentialCrossSection(const double wimpMass, const double crossSection,
                                                         const double velocity, const double recoilEnergy,
                                                         const double Anum) {
    const double formFactor = GetHelmFormFactor(recoilEnergy, Anum);

    return GetDifferentialCrossSectionNoHelmFormFactor(wimpMass, crossSection, velocity, Anum) * formFactor *
           formFactor;
}

//////////////////////////////////////////////////
/// \brief Get recoil rate for a given WIMP mass
/// and recoil energy, note that the recoil rate
/// is normalized by c/kg/day
///
const double TRestWimpUtils::GetRecoilRate(const double wimpMass, const double crossSection,
                                           const double recoilEnergy, const double Anum, const double vLab,
                                           const double vRMS, const double vEscape, const double wimpDensity,
                                           const double abundance) {
    const double vMin = GetVMin(wimpMass, Anum, recoilEnergy);
    const double vMax = vEscape + vLab;
    const double nNuclei = abundance * N_AVOGADRO * 1E3 / Anum;  // Number of atoms

    if (vMin > vMax) return 0;

    double rate = 0;
    const double velStep = 0.1;  // km/s

    for (double v = vMin; v < vMax; v += velStep) {
        const double flux = 1E5 * v * wimpDensity / wimpMass;
        const double diffRate = flux *
                                GetDifferentialCrossSection(wimpMass, crossSection, v, recoilEnergy, Anum) *
                                GetVelocityDistribution(v, vLab, vRMS, vEscape);
        rate += diffRate * velStep;
    }

    return rate * SECONDS_PER_DAY * nNuclei;
}

//////////////////////////////////////////////////
/// \brief Get Lindhard quenching factor for a
/// given recoil energy (keV) and target
///
const double TRestWimpUtils::GetQuenchingFactor(const double recoilEnergy, const double Anum,
                                                const double Znum) {
    const double deltaE = 0.0001, Emin = 0.1, resolution = 0.1;  // keV
    double g, Er = recoilEnergy, Ev;

    do {
        Er -= resolution;
        g = 0.66 * pow(Znum, 5. / 18.) / sqrt(Anum) * pow(recoilEnergy, 1. / 6.);
    } while ((Er / (1 + g) * g) > Emin);

    do {
        Er += deltaE;
        g = 0.66 * pow(Znum, 5. / 18.) / sqrt(Anum) * pow(Er, 1. / 6.);
        Ev = (Er / (1 + g) * g);
    } while (recoilEnergy > Ev);

    return (1 + g) / g;
}

//////////////////////////////////////////////////
/// \brief Parse a chemical compound into a map
/// of elements and coefficients. The compound
/// can contain parentheses to indicate a
/// subCompound. The subCompound can have a
/// coefficient after the closing parenthesis.
/// The function is recursive, so subCompounds
/// can contain subCompounds.
///
std::map<std::string, int> TRestWimpUtils::ParseChemicalCompound(const std::string& compound) {
    std::map<std::string, int> elementMap;
    std::string elementName;
    int coefficient = 1;
    
    for (size_t i = 0; i < compound.size(); ) {
        // Check for uppercase letter (start of an element)
        if (std::isupper(compound[i])) {
            elementName = compound[i];
            i++;

            // Check for lowercase letters (additional characters in element name)
            while (i < compound.size() && std::islower(compound[i])) {
                elementName += compound[i];
                i++;
            }

            // Check for a number (coefficient)
            if (i < compound.size() && std::isdigit(compound[i])) {
                coefficient = 0;
                while (i < compound.size() && std::isdigit(compound[i])) {
                    coefficient = coefficient * 10 + (compound[i] - '0');
                    i++;
                }
            }
            // Add the element and coefficient to the map
            elementMap[elementName] += coefficient;
        }
        else if (compound[i] == '(') { // Check for a subCompound inside parentheses
            i++;
            std::string subCompound;
            while (i < compound.size() && compound[i] != ')') {
                subCompound += compound[i];
                i++;
            }
            i++; // Move past the closing parenthesis
            // Find the subscript after the closing parenthesis
            coefficient = 1;
            if (i < compound.size() && std::isdigit(compound[i])) {
                coefficient = 0;
                while (i < compound.size() && std::isdigit(compound[i])) {
                    coefficient = coefficient * 10 + (compound[i] - '0');
                    i++;
                }
            }
            // Recursively call the function to handle the subCompound
            std::map<std::string, int> subElementMap = ParseChemicalCompound(subCompound);
            for (auto& pair : subElementMap) {
                elementMap[pair.first] += pair.second * coefficient;
            }
        }
        else 
            i++;
    }
    return elementMap;
}