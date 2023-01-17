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
/// functions are derived from a D. Diez-Ibañez phython code.
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

#include "TRestWimpUtils.h"

#include <TMath.h>

//////////////////////////////////////////////////
/// \brief Get relative nuclear cross section
/// within a WIMP and a nucleon, assuming
/// SCALAR INTERACTION
///
double TRestWimpUtils::GetRelativeNuclearCS(const double wimpMass, const double Anum) {
    const double reducedMass = GetReducedMass(wimpMass, Anum);
    const double reducedMassSingle = GetReducedMass(wimpMass, 1.);
    return pow(Anum * GEV_PER_UMA * reducedMass / reducedMassSingle, 2.);
}

//////////////////////////////////////////////////
/// \brief Get WIMP-nucleon reduced mass
/// (WIMP mass in GeV)
///
double TRestWimpUtils::GetReducedMass(const double wimpMass, const double Anum) {
    // WIMP mass in GeV
    return wimpMass * GEV_PER_UMA * Anum / (wimpMass + Anum * GEV_PER_UMA);
}

//////////////////////////////////////////////////
/// \brief Get Helm form factor for a given recoil
/// energy and nucleous target, mass recoil
/// energy in keV
///
double TRestWimpUtils::GetHelmFormFactor(const double recoilEnergy, const double Anum) {
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
double TRestWimpUtils::GetVMin(const double wimpMass, const double Anum, const double recoilEnergy) {
    const double reducedMass = GetReducedMass(wimpMass, Anum);
    return sqrt(Anum * GEV_PER_UMA * recoilEnergy * 1E-6 / (2. * reducedMass * reducedMass)) * LIGHT_SPEED;
}

//////////////////////////////////////////////////
/// \brief Get velocity distribution for a given
/// WIMP velocity
///
double TRestWimpUtils::GetVelocityDistribution(const double v, const double vLab, const double vRMS,
                                               const double vEscape) {
    const double vAdim = vRMS / vLab;
    const double Nesc = erf(vAdim) - (2. / sqrt(TMath::Pi())) * (vAdim)*exp(-vAdim * vAdim);
    const double xMax = std::min(1., (vEscape * vEscape - vLab * vLab - v * v) / (2. * vLab * v));
    return v * Nesc / (vLab * vRMS * sqrt(TMath::Pi())) *
           (exp(-(v - vLab) * (v - vLab) / (vRMS * vRMS)) -
            exp(-(v * v + vLab * vLab + 2 * v * vLab * xMax) / (vRMS * vRMS)));
}

//////////////////////////////////////////////////
/// \brief Differential cross section in energy,
/// in [cm/keV]. E in keV, velocity in km/s,
/// wimp mass in Gev/c^2, cross section per
/// nucleon in cm^2, Anum in atomic units (amu)
/// (or atomic weight)
///
double TRestWimpUtils::GetDifferentialCrossSection(const double wimpMass, const double crossSection,
                                                   const double velocity, const double recoilEnergy,
                                                   const double Anum) {
    const double cs = GetRelativeNuclearCS(wimpMass, Anum) * crossSection;
    const double reducedMass = GetReducedMass(wimpMass, Anum);
    const double Emax = 1E6 / LIGHT_SPEED / LIGHT_SPEED * 2. * reducedMass * reducedMass * velocity *
                        velocity / (Anum * GEV_PER_UMA);
    const double formFactor = GetHelmFormFactor(recoilEnergy, Anum);

    return cs * formFactor * formFactor / Emax;
}

//////////////////////////////////////////////////
/// \brief Get recoil rate for a given WIMP mass
/// and recoil energy, note that the recoil rate
/// is normalized by c/kg/day
///
double TRestWimpUtils::GetRecoilRate(const double wimpMass, const double crossSection,
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
double TRestWimpUtils::GetQuenchingFactor(const double recoilEnergy, const double Anum, const double Znum) {
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
