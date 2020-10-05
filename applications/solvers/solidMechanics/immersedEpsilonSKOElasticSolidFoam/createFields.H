    Info<< "Reading materials field \n" << endl;
    volScalarField materials
    (
        IOobject
        (
            "materials",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );
 
    // create a divisor that is 1e9 for all spots where materials is 1 or higher
    // and that is 1 for all spots where materials is zero (which is the fluid only)
    volScalarField materialsDivisor = (materials+scalar(1e-9))*1e9;
    //materialsDivisor.write();

    volScalarField solidDesignVariable
    (
        IOobject
        (
            "solidDesignVariable",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedScalar("zeroA", dimless, scalar(1.0))
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
    Info<< "Reading pressure field p\n" << endl;

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


    volScalarField pFluid
    (
        IOobject
        (
            "pFluid",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
	dimensionedScalar("zeroB", dimPressure, scalar(1.0))
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

    surfaceVectorField meltPressureTermf
    (
        IOobject
        (
            "meltPressureTermf",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", (dimMass*dimLength/(dimTime*dimTime)), vector::zero)
    );

    volVectorField gradMeltPressureTerm
    (
        IOobject
        (
            "gradMeltPressureTerm",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimMass/dimArea/(dimTime*dimTime), vector::zero)
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

Info<< "creating constitutive model from rheologyProperties \n" << endl;
    constitutiveModel rheology(sigma, U);
    volScalarField rho = rheology.rho();
    volScalarField mu = rheology.mu();
    volScalarField lambda = rheology.lambda();
    surfaceScalarField muf = rheology.muf();
    surfaceScalarField lambdaf = rheology.lambdaf();

Info<< "now reading information dictionary from transportProperties \n" << endl;

IOdictionary transportProperties
    (
   	IOobject
	(
		"transportProperties",
		runTime.constant(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	)
     );
    
    //lookup rho of the fluid from transport properties and make it dimensionless	
    dimensionedScalar rhoFluid(transportProperties.lookup("rhoFluid"));
    //SrhoFluid.dimensions().reset(dimless); 




    //calculate the acting fluid pressure by taking the rho-normalized pressure p, undoing the
    //normalization by multiplying with rho and take out all parts that are not located in the fluid
    //by dividing through materialsDivisor
    pFluid = p*rhoFluid/(materialsDivisor);
    //pFluid.write();

    //we look up the desired tensile strength of the materials. cells with a epsilon below our desired value will gradually set to material 2 if they were material 1 in the first place
    dimensionedScalar targetEpsilon(transportProperties.lookup("targetEpsilon"));

    dimensionedScalar timeFactor(transportProperties.lookup("timeFactor"));

    //we retain a "meltPressurePrefix" that is multiplied with pFluid before pFluid enters the matrix
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

    surfaceScalarField meltPressurePrefixf = fvc::interpolate(meltPressurePrefix, "meltPressurePrefix");
    //meltPressurePrefixf.write();

    surfaceVectorField n = mesh.Sf()/mesh.magSf();

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
