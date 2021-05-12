/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "requiredAverageCorrection.H"
#include "phasePair.H"
#include "phaseSystem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseTransferModels
{
    defineTypeNameAndDebug(requiredAverageCorrection, 0);
    addToRunTimeSelectionTable(phaseTransferModel, requiredAverageCorrection, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseTransferModels::requiredAverageCorrection::requiredAverageCorrection
(
    const dictionary& dict,
    const phasePair& pair
)
:
    phaseTransferModel(dict, pair),
    phaseName(dict.get<word>("phase")),
    reguiredAverage(dict.get<scalar>("reguiredAverage")),
    relaxationFactor(dict.getOrDefault<scalar>("relaxation",1.0)),
    totalMeshV(gSum((*(&pair_.phase1())).mesh().V())),
    cellWeigth
    (
        IOobject
        (
            "cellWeigth",
            (*(&pair_.phase1())).mesh().time().timeName(),
            (*(&pair_.phase1())).mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (*(&pair_.phase1())).mesh(),
        dimensionedScalar(dimensionSet(0,0,0,0,0,0,0), Zero),
        calculatedFvPatchField<scalar>::typeName
    )
{
	const fvMesh& mesh = (*(&pair_.phase1())).mesh();

	cellWeigth.primitiveFieldRef() = mesh.V()/totalMeshV;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseTransferModels::requiredAverageCorrection::~requiredAverageCorrection()
{}
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::tmp<Foam::volScalarField>
Foam::phaseTransferModels::requiredAverageCorrection::dmdt() const
{
    const phaseModel* phasePtr = nullptr;
    if (phaseName == pair_.first())
    {
        phasePtr = &pair_.phase1();
    }
    else if (phaseName == pair_.second())
    {
        phasePtr = &pair_.phase2();
    }
    else
    {
        FatalErrorInFunction
            << "The specified droplet phase, " << phaseName << ", is not in "
            << "the " << pair_ << " pair"
            << exit(FatalError);
    }

    const phaseModel& alpha = *phasePtr;
    const fvMesh&mesh = alpha.mesh();
    const dimensionedScalar dt(mesh.time().deltaT());
    const volScalarField& rho = alpha.rho();

    scalar currentAverage = alpha.weightedAverage(mesh.V()).value();
    scalar correction = reguiredAverage - currentAverage;
    //scalar relaxation = mesh.time().deltaTValue();
    Info << "-------------------------------------" << endl;
    Info << "currentAverage  = " << currentAverage << endl;
    Info << "requiredAverage = " << reguiredAverage << endl;
    Info << "correction      = " << correction      << endl;
    Info << "-------------------------------------" << endl;
    return relaxationFactor*cellWeigth*correction*rho/dt; //rho/dt;
}


// ************************************************************************* //
