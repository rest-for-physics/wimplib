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
/// TRestWimpNucleus is the a container class used to store different
/// nucleus parameters used in TRestWimpMetadata.
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
/// \class TRestWimpMetadata
/// \author: JuanAn Garcia   e-mail: juanangp@unizar.es
///
/// <hr>
///

#include "TRestWimpNucleus.h"
#include "TRestWimpUtils.h"

#include "TRestMetadata.h"

ClassImp(TRestWimpNucleus);

TRestWimpNucleus::TRestWimpNucleus() {
    fNucleusName = "";
    fAnum = 0;
    fZnum = 0;
    fAbundance = 0;
    fAbundanceMol = 0;
}

TRestWimpNucleus::~TRestWimpNucleus() {}

void TRestWimpNucleus::PrintNucleus() {
    RESTMetadata << "-----------------------------" << RESTendl;
    RESTMetadata << "Nuclei name " << fNucleusName.Data() << RESTendl;
    RESTMetadata << "Atomic number " << fAnum << RESTendl;
    RESTMetadata << "Number of protons " << fZnum << RESTendl;
    RESTMetadata << "Abundance " << fAbundance << RESTendl;
    RESTMetadata << "Abundance (mol) " << fAbundanceMol << RESTendl;
    RESTMetadata << "-----------------------------" << RESTendl;
}

int TRestWimpNucleus::GetStechiometricFactorFromCompound(const std::string& compound) {

    auto elementMap = TRestWimpUtils::ParseChemicalCompound(compound);

    int stechiometricFactor = 0;
    for (auto& pair : elementMap) {
        if (pair.first == fNucleusName.Data()) {
            stechiometricFactor = pair.second;
            break;
        }
    }
    if (stechiometricFactor == 0)
        RESTWarning << "No nucleus " << fNucleusName.Data() << " founnd in compound " << compound << RESTendl;
    
    return stechiometricFactor;
}