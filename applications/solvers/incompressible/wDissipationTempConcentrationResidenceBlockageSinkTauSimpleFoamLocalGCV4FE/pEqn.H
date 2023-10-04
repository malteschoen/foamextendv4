{
    p.boundaryField().updateCoeffs();

    // Prepare clean 1/Ap without contribution from under-relaxation
    // HJ, 26/Oct/2015
    volScalarField rAW
    (
        "(1|A(W))",
        1/HWEqn().A()
    );

    // Store velocity under-relaxation point before using W for
    // the flux precursor
    W.storePrevIter();

    W = rAW*HWEqn().H();
    HWEqn.clear();

    phi = fvc::interpolate(W) & mesh.Sf();
    adjustPhi(phi, W, p);

    // Non-orthogonal pressure corrector loop
    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAW, p) == fvc::div(phi)
        );

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            phi -= pEqn.flux();
        }
    }

#   include "continuityErrs.H"

    // Explicitly relax pressure for momentum corrector
    p.relax();

    // Momentum corrector
    // Note: since under-relaxation does not change aW, H/a in W can be
    // re-used.  HJ, 22/Jan/2016
    W = WUrf*(W - rAW*fvc::grad(p)) + (1 - WUrf)*W.prevIter();
    W.correctBoundaryConditions();
}
