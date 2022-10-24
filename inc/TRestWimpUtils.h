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

#include <iostream>

#ifndef RestCore_TRestWimpUtils
#define RestCore_TRestWimpUtils

/// This namespace define utilities (functions) to calculate different WIMP parameters
namespace TRestWimpUtils {

/// Physics constants
constexpr double GEV_PER_UMA = 0.93149432;
constexpr double HC_KEV_FM = 197327.053;
constexpr double LIGHT_SPEED = 300000.0;  // km/s
constexpr double SECONDS_PER_DAY = 86400;
constexpr double N_AVOGADRO = 6.0221367E23;
constexpr double MBARN_PER_GEVM2 = 0.38937966;
constexpr double CM2_PER_MBARN = 1e-27;
constexpr double FERMI_CONSTANT = 1.16639e-5;  // GeV-2

/// Generic functions for different calculations
const double GetRelativeNuclearCS(const double wimpMass, const double Anum);
const double GetReducedMass(const double wimpMass, const double Anum);
const double GetHelmFormFactor(const double recoilEnergy, const double Anum);
const double Bessel(const double x);
const double GetVMin(const double wimpMass, const double Anum, const double recoilEnergy);
const double GetVelocityDistribution(const double v, const double vLab, const double vRMS,
                                     const double vEscape);
const double GetDifferentialCrossSection(const double wimpMass, const double crossSection,
                                         const double velocity, const double recoilEnergy, const double Anum);
const double GetRecoilRate(const double wimpMass, const double crossSection, const double recoilEnergy,
                           const double Anum, const double vLab, const double vRMS, const double vEscape,
                           const double wimpDensity, const double abundance);
const double GetQuenchingFactor(const double recoilEnergy, const double Anum, const double Znum);

}  // namespace TRestWimpUtils

#endif
