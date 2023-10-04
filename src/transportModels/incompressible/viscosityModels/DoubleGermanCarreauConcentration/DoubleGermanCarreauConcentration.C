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

#include "DoubleGermanCarreauConcentration.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(DoubleGermanCarreauConcentration, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        DoubleGermanCarreauConcentration,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::DoubleGermanCarreauConcentration::calcNu() const
{
const volScalarField& Tcurrent=U_.mesh().lookupObject<volScalarField>("T");
const volScalarField& concentrationCurrent=U_.mesh().lookupObject<volScalarField>("concentration");

volScalarField dimlessConcentration = concentrationCurrent;
dimlessConcentration.dimensions().reset(dimless);

volScalarField highContribution =concentrationCurrent;
highContribution *= scalar(0.0);
highContribution.dimensions().reset(dimless);

volScalarField lowContribution =concentrationCurrent;
lowContribution *= scalar(0.0);
lowContribution.dimensions().reset(dimless);

highContribution = max(
			scalar(0),
			min(scalar(1.0), (scalar(0.5)*dimlessConcentration))
			);
lowContribution = scalar(1.0)-highContribution;

Info<< min(highContribution) << " and "<<max(highContribution)<<endl;
Info<< min(lowContribution) << " and "<<max(lowContribution)<<endl;

volScalarField nuA_local = lowContribution*nuA_for_c_zero_+highContribution*nuA_for_c_one_;
volScalarField B_local = lowContribution*B_for_c_zero_+highContribution*B_for_c_one_;
volScalarField C_local = lowContribution*C_for_c_zero_+highContribution*C_for_c_one_;
volScalarField Ts_local = lowContribution*Ts_for_c_zero_+highContribution*Ts_for_c_one_;
volScalarField Tmeasure_local = lowContribution*Tmeasure_for_c_zero_+highContribution*Tmeasure_for_c_one_;

//Info<< min(nuA_local) << " and "<<max(nuA_local)<<endl;
//Info<< min(B_local) << " and "<<max(B_local)<<endl;
//Info<< min(C_local) << " and "<<max(C_local)<<endl;
//Info<< min(Ts_local) << " and "<<max(Ts_local)<<endl;
//Info<< min(Tmeasure_local) << " and "<<max(Tmeasure_local)<<endl;

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


termB = (scalar(8.86)*(Tnew - Ts_local)) / (scalar(101.6)+(Tnew -Ts_local));
termA = (scalar(8.86)*(Tmeasure_local - Ts_local)) / (scalar(101.6)+(Tmeasure_local -Ts_local));
alpha_shift = pow(
			scalar(10),
				min(
				scalar(300),
				(termA - termB)
				)
			);


     return min
	(
			//(nuA_*alpha_shift) / pow(scalar(1)+alpha_shift*B_*strainRate(), C_) ,
			(nuA_local*alpha_shift) / pow(scalar(1)+alpha_shift*B_local*newStrainRate, C_local) ,
			(nuA_local*alpha_shift) / pow(scalar(1)+alpha_shift*B_local*newStrainRate, C_local) 


	);
	
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::DoubleGermanCarreauConcentration::DoubleGermanCarreauConcentration
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    DoubleGermanCarreauConcentrationCoeffs_(viscosityProperties.subDict(typeName + "Coeffs")),
    nuA_for_c_zero_(DoubleGermanCarreauConcentrationCoeffs_.lookup("nuA_for_c_zero")),
    nuA_for_c_one_(DoubleGermanCarreauConcentrationCoeffs_.lookup("nuA_for_c_one")),
    B_for_c_zero_(DoubleGermanCarreauConcentrationCoeffs_.lookup("B_for_c_zero")),
    B_for_c_one_(DoubleGermanCarreauConcentrationCoeffs_.lookup("B_for_c_one")),
    C_for_c_zero_(DoubleGermanCarreauConcentrationCoeffs_.lookup("C_for_c_zero")),
    C_for_c_one_(DoubleGermanCarreauConcentrationCoeffs_.lookup("C_for_c_one")),
    Ts_for_c_zero_(DoubleGermanCarreauConcentrationCoeffs_.lookup("Ts_for_c_zero")),
    Ts_for_c_one_(DoubleGermanCarreauConcentrationCoeffs_.lookup("Ts_for_c_one")),
    Tmeasure_for_c_zero_(DoubleGermanCarreauConcentrationCoeffs_.lookup("Tmeasure_for_c_zero")),
    Tmeasure_for_c_one_(DoubleGermanCarreauConcentrationCoeffs_.lookup("Tmeasure_for_c_one")),
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

bool Foam::viscosityModels::DoubleGermanCarreauConcentration::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    DoubleGermanCarreauConcentrationCoeffs_ = viscosityProperties.subDict(typeName + "Coeffs");

    DoubleGermanCarreauConcentrationCoeffs_.lookup("nuA_for_c_zero") >> nuA_for_c_zero_;
    DoubleGermanCarreauConcentrationCoeffs_.lookup("nuA_for_c_one") >> nuA_for_c_one_;

    DoubleGermanCarreauConcentrationCoeffs_.lookup("B_for_c_zero") >> B_for_c_zero_;
    DoubleGermanCarreauConcentrationCoeffs_.lookup("B_for_c_one") >> B_for_c_one_;

    DoubleGermanCarreauConcentrationCoeffs_.lookup("C_for_c_zero") >> C_for_c_zero_;
    DoubleGermanCarreauConcentrationCoeffs_.lookup("C_for_c_one") >> C_for_c_one_;

    DoubleGermanCarreauConcentrationCoeffs_.lookup("Ts_for_c_zero") >> Ts_for_c_zero_;
    DoubleGermanCarreauConcentrationCoeffs_.lookup("Ts_for_c_one") >> Ts_for_c_one_;

    DoubleGermanCarreauConcentrationCoeffs_.lookup("Tmeasure_for_c_zero") >> Tmeasure_for_c_zero_;
    DoubleGermanCarreauConcentrationCoeffs_.lookup("Tmeasure_for_c_one") >> Tmeasure_for_c_one_;


    return true;
}


// ************************************************************************* //
