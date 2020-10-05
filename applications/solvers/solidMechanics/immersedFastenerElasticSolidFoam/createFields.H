    Info<< "Reading materials field \n" << endl;
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
 
    // create a divisor that is 1e9 for all spots where materials is 1 or higher
    // and that is 1 for all spots where materials is zero (which is the fluid only)
    volScalarField materialsDivisor = (materials+scalar(1e-9))*1e9;
    //materialsDivisor.write();


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
    Info<< "Reading melt pressure field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    Info<< "Reading fastener pressure vector field pFastener\n" << endl;
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

    volTensorField gradU
    (
        IOobject
        (
            "grad(U)",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedTensor("zero", dimless, tensor::zero)
    );

    surfaceVectorField snGradU
    (
        IOobject
        (
            "snGrad(U)",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimless, vector::zero)
    );

    volVectorField V
    (
        IOobject
        (
            "V",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::ddt(U)
    );

   volTensorField gradV = fvc::ddt(gradU);
   surfaceVectorField snGradV = (snGradU - snGradU.oldTime())/runTime.deltaT();

    volSymmTensorField epsilon
    (
        IOobject
        (
            "epsilon",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    );

    volSymmTensorField sigma
    (
        IOobject
        (
            "sigma",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero)
    );

    volVectorField divSigmaExp
    (
        IOobject
        (
            "divSigmaExp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimForce/dimVolume, vector::zero)
    );

    //reading information on the mechanical properties of the solid material
    constitutiveModel rheology(sigma, U);
    volScalarField rho = rheology.rho();
    volScalarField mu = rheology.mu();
    volScalarField lambda = rheology.lambda();
    surfaceScalarField muf = rheology.muf();
    surfaceScalarField lambdaf = rheology.lambdaf();

Info<< "Reading information from transportProperties \n" << endl;

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
    
    //lookup rho of the fluid from transport properties 	
    dimensionedScalar rhoFluid(transportProperties.lookup("rhoFluid"));
    //SrhoFluid.dimensions().reset(dimless); 

    //calculate the acting fluid pressure by taking the rho-normalized pressure p, undoing the
    //normalization by multiplying with rho and take out all parts that are not located in the fluid
    //by dividing through materialsDivisor
    volScalarField pFluid = p*rhoFluid/(materialsDivisor);
    //pFluid.write();

    //we retain a "meltPressurePrefix" that is multiplied with pFluid before the results enters the matrix
    dimensionedScalar meltPressurePrefixScalar(transportProperties.lookup("meltPressurePrefixScalar"));



    volScalarField meltPressurePrefix
    (
        IOobject
        (
            "meltPressurePrefix",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
	meltPressurePrefixScalar
    );

    //meltPressurePrefix.write();
    //we construct a surfaceScalarField of the meltPressurePrefix
    surfaceScalarField meltPressurePrefixf = fvc::interpolate(meltPressurePrefix, "meltPressurePrefix");
    //meltPressurePrefixf.write();

    //we retain a "fastenerPrefix" that is multiplied with pFastener before the result enters the matrix
    dimensionedScalar fastenerPrefixScalar(transportProperties.lookup("fastenerPrefixScalar"));

    volScalarField fastenerPrefix
    (
        IOobject
        (
            "fastenerPrefix",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
	fastenerPrefixScalar
    );

    //fastenerPrefix.write();
    //we construct a surfaceScalarField of the fastenerPrefix
    surfaceScalarField fastenerPrefixF = fvc::interpolate(fastenerPrefix, "fastenerPrefix");

    surfaceVectorField n = mesh.Sf()/mesh.magSf();

    dimensionedScalar verbose(transportProperties.lookup("verbose"));


// check if there are contact boundaries to know
// if contact fields should be written for visualisation
// note: this is just for visualisation
// bool contactBoundaries = false;
// forAll(U.boundaryField(), patchi)
// {
//   if(U.boundaryField()[patchi].type()
//      ==
//      solidContactFvPatchVectorField::typeName)
//     {
//       contactBoundaries = true;
//       break;
//     }
// }


   // for aitken relaxation
   volVectorField aitkenDelta
   (
       IOobject
       (
           "aitkenDelta",
           runTime.timeName(),
           mesh,
           IOobject::NO_READ,
           IOobject::NO_WRITE
       ),
       mesh,
       dimensionedVector("zero", dimLength, vector::zero)
   );


// aitken relaxation factor
scalar aitkenInitialRes = 1.0;
scalar aitkenTheta = 0.1;
