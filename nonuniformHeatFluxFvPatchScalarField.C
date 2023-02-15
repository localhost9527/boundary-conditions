/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2020 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "nonuniformHeatFluxFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "addToRunTimeSelectionTable.H"

#include "phaseSystem.H"
#include "compressibleMomentumTransportModel.H"
#include "PhaseCompressibleMomentumTransportModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nonuniformHeatFluxFvPatchScalarField::
nonuniformHeatFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    q_(p.size(), 0.0),
    relax_(1.0),
    Tmin_(0.0)
{}


Foam::nonuniformHeatFluxFvPatchScalarField::
nonuniformHeatFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    q_("q", dict, p.size()),
    relax_(dict.lookupOrDefault<scalar>("relax", 1.0)),
    Tmin_(dict.lookupOrDefault<scalar>("Tmin", 273))
{}


Foam::nonuniformHeatFluxFvPatchScalarField::
nonuniformHeatFluxFvPatchScalarField
(
    const nonuniformHeatFluxFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(psf, p, iF, mapper),
    q_(mapper(psf.q_)),
    relax_(psf.relax_),
    Tmin_(psf.Tmin_)
{}


Foam::nonuniformHeatFluxFvPatchScalarField::
nonuniformHeatFluxFvPatchScalarField
(
    const nonuniformHeatFluxFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(psf, iF),
    q_(psf.q_),
    relax_(psf.relax_),
    Tmin_(psf.Tmin_)
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::nonuniformHeatFluxFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Lookup the fluid model
    const phaseSystem& fluid =
        (
            db().lookupObject<phaseSystem>("phaseProperties")
        );

    const scalarField& Tp = *this;
    //scalarField& Tp = *this;

    scalarField A(Tp.size(), scalar(0));
    scalarField B(Tp.size(), scalar(0));
    scalarField Q(Tp.size(), scalar(0));

    forAll(fluid.phases(), phasei)
    {
        const phaseModel& phase = fluid.phases()[phasei];
        const fluidThermo& thermo = phase.thermo();

        const fvPatchScalarField& alpha =
            phase.boundaryField()[patch().index()];

        const fvPatchScalarField& T =
            thermo.T().boundaryField()[patch().index()];

        const scalarField kappaEff(phase.kappaEff(patch().index()));

        if (debug)
        {
            scalarField q0(T.snGrad()*alpha*kappaEff);
            Q += q0;

            Info<< patch().name() << " " << phase.name()
                << ": Heat flux " << gMin(q0) << " - " << gMax(q0) << endl;
        }

        A += T.patchInternalField()*alpha*kappaEff*patch().deltaCoeffs();
        B += alpha*kappaEff*patch().deltaCoeffs();
    }

    if (debug)
    {
        Info<< patch().name() << " " << ": overall heat flux "
            << gMin(Q) << " - " << gMax(Q) << " W/m2, power: "
            << gSum(patch().magSf()*Q) << " W" << endl;
    }
    //user defined
    const vectorField& Cf = patch().Cf();   //center of the patch
    scalarField Shapefunc(Tp.size(), scalar(0));
    constexpr scalar 	pi (M_PI);
    forAll(Cf, faceI)
    {
    	const scalar x = Cf[faceI][0];
    	//operator==((1 - relax_)*Tp[faceI] + relax_*max(Tmin_,(q_ + A)/(B)));
    	Shapefunc[faceI]=Foam::sin(x/2*pi);
    }
    operator==((1 - relax_)*Tp + relax_*max(Tmin_,(q_*Shapefunc + A)/(B)));

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::nonuniformHeatFluxFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    m(q_, q_);
}


void Foam::nonuniformHeatFluxFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const nonuniformHeatFluxFvPatchScalarField& mptf =
        refCast<const nonuniformHeatFluxFvPatchScalarField>(ptf);

    q_.rmap(mptf.q_, addr);
}


void Foam::nonuniformHeatFluxFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeEntry(os, "relax", relax_);
    writeEntry(os, "q", q_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        nonuniformHeatFluxFvPatchScalarField
    );
}


// ************************************************************************* //
