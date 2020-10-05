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


volVectorField U
(
    IOobject
    (
        "U",
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

volTensorField outputGradU
(
    IOobject
    (
        "outputGradU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedTensor("outputGradU", dimensionSet(0,0,-1,0,0,0,0), tensor::zero)
);


volScalarField detgradU
(
    IOobject
    (
        "detgradU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("detgradU", dimensionSet(0,0,0,0,0,0,0), scalar (0.0))
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



