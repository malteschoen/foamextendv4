/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "GermanCarreau.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(GermanCarreau, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        GermanCarreau,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::GermanCarreau::calcNu() const
{
const volScalarField& Tcurrent=U_.mesh().lookupObject<volScalarField>("T");

volScalarField alpha_shift =Tcurrent;			//instead of constructing a new volScalarField, we "copy" the old one			
alpha_shift *= scalar(0.0);				//we then set all values of the volScalarField to zero
alpha_shift.dimensions().reset(dimless);		//also, since T is a Value in Kelvin, we make the new Field dimensionsless

volScalarField termA =Tcurrent;
termA *= scalar(0.0);
termA.dimensions().reset(dimless);

volScalarField termB =Tcurrent;
termB *= scalar(0.0);
termB.dimensions().reset(dimless);

volScalarField Tnew =Tcurrent;
Tnew.dimensions().reset(dimless);

dimensionSet arDims (0,0,-1,0,0,0,0);
volScalarField newStrainRate = Tcurrent;
newStrainRate.dimensions().reset(arDims);
newStrainRate *= scalar(0.0);
newStrainRate = strainRate();
newStrainRate.dimensions().reset(dimless);
newStrainRate = newStrainRate + SMALL;
newStrainRate.dimensions().reset(arDims);


termB = (scalar(8.86)*(Tnew - Ts_)) / (scalar(101.6)+(Tnew -Ts_));
termA = (scalar(8.86)*(Tmeasure_ - Ts_)) / (scalar(101.6)+(Tmeasure_ -Ts_));
alpha_shift = pow(scalar(10), (termA - termB));


     return min
	(
			//(nuA_*alpha_shift) / pow(scalar(1)+alpha_shift*B_*strainRate(), C_) ,
			(nuA_*alpha_shift) / pow(scalar(1)+alpha_shift*B_*newStrainRate, C_) ,
			(nuA_*alpha_shift) / pow(scalar(1)+alpha_shift*B_*newStrainRate, C_) 


	);
	
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::GermanCarreau::GermanCarreau
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    GermanCarreauCoeffs_(viscosityProperties.subDict(typeName + "Coeffs")),
    nuA_(GermanCarreauCoeffs_.lookup("nuA")),
    B_(GermanCarreauCoeffs_.lookup("B")),
    C_(GermanCarreauCoeffs_.lookup("C")),
    Ts_(GermanCarreauCoeffs_.lookup("Ts")),
    Tmeasure_(GermanCarreauCoeffs_.lookup("Tmeasure")),
    nu_
	(
        IOobject
		(
			name,
			U_.time().timeName(),
			U_.db(),
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		
		),	
       		calcNu()
	)
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::GermanCarreau::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    GermanCarreauCoeffs_ = viscosityProperties.subDict(typeName + "Coeffs");

    GermanCarreauCoeffs_.lookup("nuA") >> nuA_;
    GermanCarreauCoeffs_.lookup("B") >> B_;
    GermanCarreauCoeffs_.lookup("C") >> C_;
    GermanCarreauCoeffs_.lookup("Ts") >> Ts_;
    GermanCarreauCoeffs_.lookup("Tmeasure") >> Tmeasure_;



    return true;
}


// ************************************************************************* //
