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

#include "localGermanCarreau.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(localGermanCarreau, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        localGermanCarreau,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::localGermanCarreau::calcNu() const
{
const volScalarField& Tcurrent=U_.mesh().lookupObject<volScalarField>("T");
const volScalarField& localnuA=U_.mesh().lookupObject<volScalarField>("localnuA");
const volScalarField& localB=U_.mesh().lookupObject<volScalarField>("localB");
const volScalarField& localC=U_.mesh().lookupObject<volScalarField>("localC");
const volScalarField& localTs=U_.mesh().lookupObject<volScalarField>("localTs");
const volScalarField& localTmeasure=U_.mesh().lookupObject<volScalarField>("localTmeasure");

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


termB = (scalar(8.86)*(Tnew - localTs)) / (scalar(101.6)+(Tnew -localTs));
termA = (scalar(8.86)*(localTmeasure - localTs)) / (scalar(101.6)+(localTmeasure -localTs));
alpha_shift = pow(
			scalar(10),
				min(
				scalar(300),
				(termA- termB)
				)
			);


     return min
	(

			(localnuA*alpha_shift) / pow(scalar(1)+alpha_shift*localB*newStrainRate, localC) ,
			(localnuA*alpha_shift) / pow(scalar(1)+alpha_shift*localB*newStrainRate, localC)


	);
	
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::localGermanCarreau::localGermanCarreau
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    localGermanCarreauCoeffs_(viscosityProperties.subDict(typeName + "Coeffs")),
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

bool Foam::viscosityModels::localGermanCarreau::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    localGermanCarreauCoeffs_ = viscosityProperties.subDict(typeName + "Coeffs");

    return true;
}


// ************************************************************************* //
