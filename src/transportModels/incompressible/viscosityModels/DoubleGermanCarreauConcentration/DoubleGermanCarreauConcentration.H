/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Class
    Foam::viscosityModels::DoubleGermanCarreauConcentration

Description
    An incompressible isothermal Carreau Three-factor viscosity model with WLF correction for temperatures

SourceFiles
    DoubleGermanCarreauConcentration.C

\*---------------------------------------------------------------------------*/

#ifndef DoubleGermanCarreauConcentration_H
#define DoubleGermanCarreauConcentration_H

#include "viscosityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{

/*---------------------------------------------------------------------------*\
                  Class DoubleGermanCarreauConcentration Declaration
\*---------------------------------------------------------------------------*/

class DoubleGermanCarreauConcentration
:
    public viscosityModel
{
    // Private data

        dictionary DoubleGermanCarreauConcentrationCoeffs_;

        dimensionedScalar nuA_for_c_zero_;
        dimensionedScalar nuA_for_c_one_;
        dimensionedScalar B_for_c_zero_;
        dimensionedScalar B_for_c_one_;
        dimensionedScalar C_for_c_zero_;
        dimensionedScalar C_for_c_one_;
        dimensionedScalar Ts_for_c_zero_;
        dimensionedScalar Ts_for_c_one_;
        dimensionedScalar Tmeasure_for_c_zero_;
        dimensionedScalar Tmeasure_for_c_one_;


        volScalarField nu_;



    // Private Member Functions

        //- Calculate and return the laminar viscosity
        tmp<volScalarField> calcNu() const;


public:

    //- Runtime type information
    TypeName("DoubleGermanCarreauConcentration");


    // Constructors

        //- Construct from components
        DoubleGermanCarreauConcentration
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    ~DoubleGermanCarreauConcentration()
    {}


    // Member Functions

        //- Return the laminar viscosity
        virtual const volScalarField& nu() const
        {
            return nu_;
        }

        //- Correct the laminar viscosity
        void correct()
        {
            nu_ = calcNu();
        }

        //- Read transportProperties dictionary
        bool read(const dictionary& viscosityProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
