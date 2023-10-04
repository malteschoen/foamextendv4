{    

   dimensionedScalar source
    (
            "source", dimensionSet(0,0,0,0,0,0,0),1.0
    );

    fvScalarMatrix residenceTimeEqn
    (
	fvm::div(phi, residenceTime)

      == 
	fvm::laplacian(DresidenceTime, residenceTime)
	+ source
    );

    residenceTimeEqn.relax();
    residenceTimeEqn.solve().initialResidual();
}



