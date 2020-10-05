Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

dimensionedScalar rho
(
    "rho",
    dimDensity,
    transportProperties.lookup("rho")
);

Info<< "Reading fields\n" << endl;

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


volScalarField nu
(
    IOobject
    (
        "nu",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

/*IOobject nuHeader
    (
        "nu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );

volScalarField* fieldNuPointer;
bool fieldNuExist = false;

if (nuHeader.headerOk())
{
fieldNuExist = true;
Info<< "Reading field nu\n" << endl;

fieldNuPointer = new volScalarField
(
    IOobject
    (
        "nu",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
);
}
*/


volVectorField W
(
    IOobject
    (
        "W",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField p_bar
(
    IOobject
    (
        "p_bar",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("p_bar", dimPressure, scalar (0.0))
);

volVectorField tauFixed
(
    IOobject
    (
        "tauFixed",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("tauFixed", dimPressure, vector::zero)
);

volScalarField strainRate
(
    IOobject
    (
        "strainRate",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("strainRate", dimensionSet(0,0,-1,0,0,0,0), scalar (0.0))
);


volVectorField elongRates
(
    IOobject
    (
        "elongRates",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("elongRates", dimensionSet(0,0,-1,0,0,0,0), vector::zero)
);

volTensorField outputGradW
(
    IOobject
    (
        "outputGradW",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedTensor("outputGradW", dimensionSet(0,0,-1,0,0,0,0), tensor::zero)
);



volScalarField eta
(
    IOobject
    (
        "eta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("eta", dimensionSet(1,-1,-1,0,0,0,0), scalar (0.0))
);



