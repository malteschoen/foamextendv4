epsilon = symm(gradU);

vector Iv(1,1,1);

symmTensor EI(1,0,0,1,0,1);

//symmTensor pStuff(pFastener.component(0),0,0,pFastener.component(1),0,pFastener.component(2));

sigma = 2*mu*epsilon + lambda*(I*tr(epsilon)) -I*meltPressurePrefix*pFluid - cmptMultiply(symm(pFastener*Iv),EI)*fastenerPrefix;

//sigma = 2*mu*epsilon + lambda*(I*tr(epsilon)) -I*meltPressurePrefix*pFluid - pStuff*fastenerPrefix;

