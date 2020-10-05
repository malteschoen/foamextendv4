/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

Application
    extrusionData

Description
    Converts the rho-normalized pressure p to a regular Pressure in bars.
    Calculates shear rates.
    Converts kinematic viscosity into dynamic viscosity.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"


    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



//vectorizedTau = *tau;
vector Iv(1,1,1);

//unityVector.write();
//volScalarField materialsDivisor = (materials+scalar(1e-9))*1e9;
vectorizedTau = (tau&unityVector)*materials;
vectorizedTau.write();

volVectorField tauStar ("tauStar", vectorizedTau);
tauStar = *tau;
tauStar.write();

volVectorField tauDiv ("tauDiv", vectorizedTau);
dimensionSet arDims1 (0,1,-2,0,0,0,0);
tauDiv.dimensions().reset(arDims1);
tauDiv = fvc::div(tau);
tauDiv.write();


dimensionSet dimLessSet (0,0,0,0,0,0,0);

meshVolumes.internalField() = mesh.V();
meshVolumes.dimensions().reset(dimLessSet);
meshVolumes.write();

volScalarField meshLengths ("meshLengths",meshVolumes);

meshLengths.dimensions().reset(dimLessSet);
meshLengths = pow(meshVolumes, 0.33333);
dimensionSet arDimsMLX (0,0,1,0,0,0,0);
meshLengths.dimensions().reset(arDimsMLX);
meshLengths.write();

volVectorField tauDivLength("tauDivLength", tauDiv);
dimensionSet arDims99 (0,1,-1,0,0,0,0);
tauDivLength.dimensions().reset(arDims99);
tauDivLength = tauDiv*meshLengths;
tauDivLength.write();

volTensorField gradW = fvc::grad(W);
volTensorField gradWT = gradW.T();


divGradW = fvc::div(gradW);
divGradW.write();

divSymmGradW = fvc::div(symm(gradW));
divSymmGradW.write();

normalVector = mesh.Sf()/mesh.magSf();
normalVector.write();


dimensionSet arDims2 (0,2,-2,0,0,0,0);
dimensionSet arDims3 (0,0,-1,0,0,0,0);

volVectorField tauDirect ("tauDirect", vectorizedTau);
	tauDirect.dimensions().reset(arDims2);
    	//tauDirect = scalar(2.0)*(*(symm(gradWT)));
	tauDirect = nu*(Iv & symm(gradW));
    	tauDirect.write();



tauF = fvc::interpolate(tau);
tauF.write();

dimensionSet arDims4 (0,0,0,0,0,0,0);

surfaceScalarField tempSSF = scalar(1e16)*phi;
   
//Info<< "break one\n" << endl;

tempSSF.dimensions().reset(arDims4);
tempSSF.write();

//Info<< "break two\n" << endl;

tannedPhi.dimensions().reset(arDims4);
tannedPhi = tanh(tempSSF);
tannedPhi.write();

//surfaceVectorField gradTannedPhi = fvc::grad(tannedPhi);

volVectorField tempW ("tempW", W);

tempW = W*scalar(1e12);
tempW.dimensions().reset(arDims4);
tempW.write();

/*tannedW = vector(
		tanh(tempW.component(vector::X)),
		tanh(tempW.component(vector::Y)),
		tanh(tempW.component(vector::Z))
		);
tannedW.write();
*/

tannedW.internalField().replace(vector::X, tanh(tempW.component(vector::X)));
tannedW.internalField().replace(vector::Y, tanh(tempW.component(vector::Y)));
tannedW.internalField().replace(vector::Z, tanh(tempW.component(vector::Z)));
tannedW.write();

volVectorField magTauTannedW ("magTauTannedW", vectorizedTau);
dimensionSet arDims5 (0,2,-2,0,0,0,0);
magTauTannedW.dimensions().reset(arDims5);
magTauTannedW = mag(tau)*tannedW;
magTauTannedW.write();
	
tauLikePhi =  tauF & normalVector;
//tauLikePhi =   normalVector & tauF ;
tauLikePhi.write();



//
//surfaceVectorField fuckyTerm = rho*(fvc::interpolate(vectorizedTau)*mesh.magSf());
//fuckyTerm.write();

//tauForceTerm = rho*cmptMultiply(
//				fvc::interpolate(vectorizedTau),mesh.Sf()
//
//				);
//
//tauForceTerm = 	(
//				fvc::interpolate(tauDiv)
//				*mesh.magSf()				
//				//^
//				//mesh.Sf()
//			);


//tauForceTerm.write();

//nuRhoWForceTerm = fvc::snGrad(W);//*rho*nu;
//nuRhoWForceTerm.write();





						

   Info<< "Action complete\n" << endl;

    return 0;
}


// ************************************************************************* //
