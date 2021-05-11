/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::phaseTransferModels::requiredAverageCorrection

Description
    Phase transfer model representing change from a dispersed phase to a film as
    a result of requiredAverageCorrection onto a third phase

SourceFiles
    requiredAverageCorrection.C

\*---------------------------------------------------------------------------*/

#ifndef requiredAverageCorrection_H
#define requiredAverageCorrection_H

#include "phaseTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace phaseTransferModels
{

/*---------------------------------------------------------------------------*\
                          Class requiredAverageCorrection Declaration
\*---------------------------------------------------------------------------*/

class requiredAverageCorrection
:
    public phaseTransferModel
{
private:

    // Private data

        //- The name of the phase which deposits
        const word phaseName;

        //- The requiredAverageCorrection efficiency
        const scalar reguiredAverage;

        const scalar relaxationFactor;

        const scalar totalMeshV;

public:

    //- Runtime type information
    TypeName("requiredAverageCorrection");


    // Constructors

        //- Construct from components
        requiredAverageCorrection
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~requiredAverageCorrection();


    // Member Functions

        //- The mass transfer rate
        virtual tmp<volScalarField> dmdt() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace phaseTransferModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
