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

    Info<<"now updating materials according to SKO principle" << endl;

    localWork =
mag(sigma.component(tensor::XX)*U.component(vector::X))+
mag(sigma.component(tensor::YY)*U.component(vector::Y))+
mag(sigma.component(tensor::ZZ)*U.component(vector::Z));

    //we slowly increase the target epsilon by multiplying with numbers inching closer to 1 with every time step

     dimensionedScalar newTargetWork = targetWork*
	min(scalar(1.0), runTime.value()/timeFactor.value());

    Info<< newTargetWork << endl;		
     int nrChanges =0;

     forAll(localWork, iterator)
     {
     nrChanges =0;
     	if (localWork[iterator] < newTargetWork.value() && changeAllowed[iterator] == 1)
		{
		materials[iterator] = scalar(2.0);
                nrChanges = nrChanges +1;
		}
     }
    Info<< nrChanges << endl;		

    runTime.write();
}
