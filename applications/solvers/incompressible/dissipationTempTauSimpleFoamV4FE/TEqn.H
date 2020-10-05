{    

    //Credit goes to Martin Becker of DHCAE Tools (martin_becker@dhcae-tools.de)	 

    volTensorField gradU = fvc::grad(U);
    volScalarField nu = laminarTransport.nu();
    tau = nu * (gradU + gradU.T()); 

    //volScalarField increaser = laminarTransport.nu();

   // volScalarField DTeffective = scalar(1.0)*alpha;

    fvScalarMatrix TEqn
    (
	fvm::div(phi, T)
	- fvm::laplacian(DT, T)

      == 
	 (1/c)*(tau && gradU)
       
    );

    TEqn.relax();
    TEqn.solve().initialResidual();

}



