    Info<< "\nReading transportProperties\n" << endl;

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

    Info<< "Reading scalar for thermal diffusivity alpha_t\n" << endl;
    dimensionedScalar DT
    (
        transportProperties.lookup("DT")
    );


    Info<< "Reading scalar for specific heat capacity c\n" << endl;
    dimensionedScalar c 
    (
        transportProperties.lookup("c")
    );

    dimensionedScalar rhoFluid
    (
        transportProperties.lookup("rho")
    );


    Info<< "Reading field T\n" << endl;
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


    Info<< "Reading field p\n" << endl;
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


    Info<< "Reading field U\n" << endl;
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

#   include "createPhi.H"

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, piso.dict(), pRefCell, pRefValue);
    mesh.schemesDict().setFluxRequired(p.name());

    singlePhaseTransportModel fluid(U, phi);

//singlePhaseTransportModel laminarTransport(U,phi);
//autoPtr<incompressible::turbulenceModel> turbulence
//(
//incompressible::turbulenceModel::New(U,phi,laminarTransport)
//);

//     pRefCell = mesh.findNearestCell(vector(0,0.205,0));
//     Info << "Ref. cell: " << pRefCell << endl;
//     Info << "Ref. cell position: " << mesh.C()[pRefCell] << endl;
