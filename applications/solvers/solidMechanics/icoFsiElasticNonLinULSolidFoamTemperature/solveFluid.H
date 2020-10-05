
        if(mesh.moving())
        {
            // Make the fluxes relative
            phi -= fvc::meshPhi(U);
        }

#       include "CourantNo.H"
#	include "singlePhaseTransportModel.H"



	/*
	const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
	(
		IOobject::groupName
		(
			turbulenceModel::propertiesName,
			internalField().group()
		)
	);

	volScalarField nu  = turbModel.nu()
	
	*/

	volScalarField nu  = fluid.nu();

        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          - fvm::laplacian(nu, U)
        );

        solve(UEqn == -fvc::grad(p));

        // --- PISO loop
        volScalarField rUA = 1.0/UEqn.A();

        while (piso.correct())
        {
            U = rUA*UEqn.H();
            phi = (fvc::interpolate(U) & mesh.Sf());

            adjustPhi(phi, U, p);

            while (piso.correctNonOrthogonal())
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rUA, p)
                 == fvc::div(phi)
                );

                pEqn.setReference(pRefCell, pRefValue);
                pEqn.solve();

                if (piso.finalNonOrthogonalIter())
                {
                    phi -= pEqn.flux();
                }
            }

#           include "continuityErrs.H"

            U -= rUA*fvc::grad(p);
            U.correctBoundaryConditions();

	   }
	 volTensorField gradU = fvc::grad(U);
	 volTensorField tau = nu * (gradU + gradU.T());

	fvScalarMatrix TEqn
	(
	fvm::ddt(T) - fvm::div(phi, T)
	==
	fvm::laplacian(DT, T) + (1/c)*(tau && gradU)
	);
	
