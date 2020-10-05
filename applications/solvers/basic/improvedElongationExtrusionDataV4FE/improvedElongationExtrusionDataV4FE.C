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



p_bar = p*(1e-5)*rho;
p_bar.write();
    Info<< "\nCalculated p_bar" << endl;

strainRate = scalar(1.41421)*mag(symm(fvc::grad(U)));
strainRate.write();
    Info<< "\nCalculated shearRate" << endl;

outputGradU = fvc::grad(U);
outputGradU.write();
    Info<< "\nCalculated gradU" << endl;

eta = nu*rho;
eta.write();	
    Info<< "\nCalculated dynamic viscosities eta" << endl;


tmp<volTensorField> gradUcalledD = fvc::grad(U);
    Info<< "\nloaded temp tensor field spoint 1" << endl;

volTensorField dimlessOutputGradU = outputGradU;
dimlessOutputGradU.dimensions().reset(dimless);
detgradU = det(dimlessOutputGradU);
detgradU.write();
    Info<< "\nCalculated dimensionless det(GradU) (positive parts only)" << endl;

						

   Info<< "Calculation of elongational extrusion data complete\n" << endl;

    return 0;
}


// ************************************************************************* //
