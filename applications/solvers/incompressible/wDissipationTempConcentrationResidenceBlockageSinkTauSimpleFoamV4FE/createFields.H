Info<< "Reading materials\n" << endl;
volScalarField materials
(
    IOobject
    (
        "materials",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading thermal diffusivity field DT\n" << endl;
volScalarField DT
(
    IOobject
    (
        "DT",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading residenceTime\n" << endl;
volScalarField residenceTime
(
    IOobject
    (
        "residenceTime",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading rho-normalized pressure field p\n" << endl;
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


Info<< "Reading velocity field W\n" << endl;
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

Info<< "Reading temperature field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading deflection field U\n" << endl;
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

Info<< "Reading concentration field\n" << endl;
volScalarField concentration
(
    IOobject
    (
        "concentration",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading sink field\n" << endl;
volScalarField sink
(
    IOobject
    (
        "sink",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volTensorField tau
(
    IOobject
    (
        "tau",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedTensor("tau", dimensionSet(0,2,-2,0,0,0,0), tensor::zero)
);

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


Info<< "Reading scalar for fluid specific heat capacity c\n" <<endl;
    dimensionedScalar c
    (
        transportProperties.lookup("c")
    );

Info<< "Reading scalar DresidenceTime for diffusion of residenceTime to increase stability\n" <<endl;
    dimensionedScalar DresidenceTime
    (
        transportProperties.lookup("DresidenceTime")
    );


Info<< "Reading scalar materialsToAlpha which converts materials to alpha\n" <<endl;
    dimensionedScalar materialsToAlpha
    (
        transportProperties.lookup("materialsToAlpha")
    );

Info<< "Reading scalar Dconcentration for concentration diffusivity\n" <<endl;
    dimensionedScalar Dconcentration
    (
        transportProperties.lookup("Dconcentration")
    );


Info<< "Reading/calculating face flux field phi\n" << endl;
surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(W) & mesh.Sf()
);



label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, simple.dict(), pRefCell, pRefValue);
mesh.schemesDict().setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(W, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(W, phi, laminarTransport)
);

