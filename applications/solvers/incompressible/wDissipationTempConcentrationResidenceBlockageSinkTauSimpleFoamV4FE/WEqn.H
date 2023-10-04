    // Solve the momentum equation

    tmp<fvVectorMatrix> HWEqn
    (
        fvm::div(phi, W)
      + turbulence->divDevReff()
      + fvm::Sp((materials*materialsToAlpha),W)
    );

    // Get under-relaxation factor
    const scalar WUrf = mesh.solutionDict().equationRelaxationFactor(W.name());

    // Momentum solution
    solve
    (
        relax(HWEqn(), WUrf)
     ==
        -fvc::grad(p)
    );
