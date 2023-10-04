{    

    fvScalarMatrix concentrationEqn
    (
	fvm::div(phi, concentration)

      == 
	fvm::laplacian(Dconcentration, concentration)
    );

    concentrationEqn.relax();
    concentrationEqn.solve().initialResidual();
}



