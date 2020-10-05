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
    transportProperties.lookup("rho")
);





Info<< "Reading fields\n" << endl;

volTensorField tau
(
    IOobject
    (
        "tau",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

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

volVectorField pFastener
(
    IOobject
    (
        "pFastener",
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

volScalarField meshVolumes
(
    IOobject
    (
        "meshVolumes",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0, 3, 0, 0, 0, 0, 0), scalar(0))
);

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

volVectorField unityVector
(
    IOobject
    (
        "unityVector",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("unity", dimless, vector(1,1,1))
);
volVectorField vectorizedTau
(
    IOobject
    (
        "vectorizedTau",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimensionSet(0, 2, -2, 0, 0, 0, 0), vector::zero)
);

volVectorField divGradW
(
    IOobject
    (
        "divGradW",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimensionSet(0, -1, -1, 0, 0, 0, 0), vector::zero)
);

volVectorField divSymmGradW
(
    IOobject
    (
        "divSymmGradW",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimensionSet(0, -1, -1, 0, 0, 0, 0), vector::zero)
);



surfaceVectorField tauForceTerm
(
    IOobject
    (
        "tauForceTerm",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimensionSet(1, 0, -2, 0, 0, 0, 0), vector::zero)

);

surfaceVectorField tauLikePhi
(
    IOobject
    (
        "tauLikePhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimensionSet(0, 2, -2, 0, 0, 0, 0), vector::zero)

);

surfaceScalarField tannedPhi
(
    IOobject
    (
        "tannedPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimensionSet(0, 3, -1, 0, 0, 0, 0), scalar(0))

);

volVectorField tannedW
(
    IOobject
    (
        "tannedW",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimensionSet(0, 0, 0, 0, 0, 0, 0), vector::zero)
);


surfaceVectorField nuRhoWForceTerm
(
    IOobject
    (
        "nuRhoWForceTerm",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimensionSet(0, 0, -1, 0, 0, 0, 0), vector::zero)

);

surfaceVectorField normalVector
(
    IOobject
    (
        "normalVector",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimensionSet(0, 0, 0, 0, 0, 0, 0), vector::zero)

);

surfaceVectorField surfaceAreaNormalVector
(
    IOobject
    (
        "surfaceAreaNormalVector",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimensionSet(0, 0, 0, 0, 0, 0, 0), vector::zero)

);

surfaceTensorField tauF
(
    IOobject
    (
        "tauF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedTensor("zero", dimensionSet(0, 2, -2, 0, 0, 0, 0), tensor::zero)

);

surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
mesh
);





