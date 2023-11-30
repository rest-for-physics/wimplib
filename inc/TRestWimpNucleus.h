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

#ifndef RestCore_TRestWimpNucleus
#define RestCore_TRestWimpNucleus

#include <TString.h>

#include <iostream>

/// A class to store different nucleus parameters
class TRestWimpNucleus {
   public:
    /// Nucleus name
    TString fNucleusName;
    /// Atomic number in amus
    Double_t fAnum;
    /// Number of protons
    Int_t fZnum;
    /// Abundance, in mass percentage
    Double_t fAbundance;
    /// Abundance, in mole (or volume)
    Double_t fAbundanceMol;

    void PrintNucleus();
    int GetStechiometricFactorFromCompound(const std::string& compound);

    // Constructor
    TRestWimpNucleus();

    // Destructor
    virtual ~TRestWimpNucleus();

    ClassDef(TRestWimpNucleus, 2);
};

#endif
