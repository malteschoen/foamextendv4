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

    Info<<"now removing material according to SKO principle" << endl;

    dimensionedScalar newTargetSigma = targetSigma*min(
                                  runTime.value()/timeFactor.value(),
                                   scalar(1.0)
                                   );
    Info<< newTargetSigma << endl;

     materials = max (scalar(0),
		materials-
		min(
		    scalar(1.0),
		    pos((newTargetSigma/sigmaEq)-scalar(2.0))
                    )
		);
    Info<<"writing down" << endl;
    runTime.write();
}
