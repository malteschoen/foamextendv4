if (runTime.outputTime())
{
    volScalarField epsilonEq
    (
        IOobject
        (
            "epsilonEq",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sqrt((2.0/3.0)*magSqr(dev(epsilon)))
    );

    Info<< "Max epsilonEq = " << max(epsilonEq).value()
        << endl;

    volScalarField sigmaEq
    (
        IOobject
        (
            "sigmaEq",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sqrt((3.0/2.0)*magSqr(dev(sigma)))
    );

    Info<< "Max sigmaEq = " << max(sigmaEq).value()
        << endl;

    Info<<"now updating solidDesignVariable according to SKO principle" << endl;

    //we slowly increase the targetSigma by multiplying with numbers inching closer to 1 with every time step
    dimensionedScalar newTargetSigma = targetSigma*min(
                                  runTime.value()/timeFactor.value(),
                                   scalar(1.0)
                                   );
    Info<< newTargetSigma << endl;


     solidDesignVariable = newTargetSigma/sigmaEq;
			

     forAll(solidDesignVariable, iterator)
     {
     	if (solidDesignVariable[iterator] < scalar(0.5))
		{
		materials[iterator] = scalar(2.0);
		}
     }

     //Smaterials = materials- pos(solidDesignVariable;

    runTime.write();
}
