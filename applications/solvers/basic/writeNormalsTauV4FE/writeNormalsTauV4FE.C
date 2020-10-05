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

//vector Iv(1,1,1);

dimensionSet dimAccelSet (0,1,-2,0,0,0,0);
dimensionSet dimLessSet (0,0,0,0,0,0,0);
dimensionSet dimMeterSet (0,0,1,0,0,0,0);
dimensionSet arDims99 (1,-2,-1,0,0,0,0);



//calculation of the divergence of tau
volVectorField tauDiv ("tauDiv", vectorizedTau);
tauDiv.dimensions().reset(dimAccelSet);
tauDiv = fvc::div(tau);
//tauDiv.write();

//calculation of cell volumina 
meshVolumes.internalField() = mesh.V();
meshVolumes.dimensions().reset(dimLessSet);
//meshVolumes.write();


//calculation of characteristical cell length
volScalarField meshLengths ("meshLengths",meshVolumes);
meshLengths.dimensions().reset(dimLessSet);
meshLengths = pow(meshVolumes, 0.33333);
meshLengths.dimensions().reset(dimMeterSet);
meshLengths.write();

//calculation of shear stresses from divergence of tau
//by multiplying it with the characteristical cell length
volVectorField tauDivLength("tauDivLength", tauDiv);
tauDivLength.dimensions().reset(arDims99);
tauDivLength = tauDiv*meshLengths*rho;
tauDivLength.write();

//calculation of normals
normalVector = mesh.Sf()/mesh.magSf();
normalVector.write();

//calculation of an alternative velocity that only indicates velocity direction by values of -1/+1
volVectorField tempW ("tempW", W);
tempW = W*scalar(1e12);
tempW.dimensions().reset(dimLessSet);

volVectorField directionW ("directionW", W);
directionW.dimensions().reset(dimLessSet);
directionW.internalField().replace(vector::X, tanh(tempW.component(vector::X)));
directionW.internalField().replace(vector::Y, tanh(tempW.component(vector::Y)));
directionW.internalField().replace(vector::Z, tanh(tempW.component(vector::Z)));
directionW.dimensions().reset(dimVelocity);					

   Info<< "Action complete\n" << endl;

    return 0;
}


// ************************************************************************* //
