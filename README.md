# foamextendv4
this was built with foam-extend version 4.0. however, it also runs under fe4.1 and fe5.0

# what is this?
we have solvers for extrusion operation (stationary, incompressible) taking into account both solids and fluids. More precisely, we have ...
1. ... immersed body solid mechanics taking into account the fluid forces (pressure, shear stress). This is a two-step, non-coupled approach - i.e. after a flow calculation there is no feedback of the geometry deformation into the fluid simulation. 
2. ... immersed body fluid flow in which fluid and solid domain share a heat/temperature equation. The funky stuff is the addition of heat sinks.
3. ... a coextrusion model based on concentration tracking - depending on concentration, different values of the Carreau-WLF-model are used in each cell

# i just wanna try solid mechanics
* splendid!
* before you go ahead, you'll need to compile wDissipationTempConcentrationResidenceBlockageSinkTauSimpleFoamV4FE
  * go to 'applications/solvers/incompressible/wDissipationTempConcentrationResidenceBlockageSinkTauSimpleFoamV4FE'
  * run 'wmake'
* also, you will need to compile immersedTauFastenerElasticSolidFoam 
  * go to '/applications/solvers/solidMechanics/immersedTauFastenerElasticSolidFoam'
  * run 'wmake'
* then head for the tutorialSolidMechanics folder
* having arrived there, follow the local readMe instructions


# i just wanna try the heat-sink stuff
* excellent!
* before you go ahead, you'll need to compile wDissipationTempConcentrationResidenceBlockageSinkTauSimpleFoamV4FE.
  *  go to 'applications/solvers/incompressible/wDissipationTempConcentrationResidenceBlockageSinkTauSimpleFoamV4FE'
  *  run 'wmake'
* then head for the tutorialHeatSink folder
* having arrived there, follow the local readMe instructions

# i just wanna try coextrusion
* magnificient!
* before you go ahead, you'll need to compile wDissipationTempConcentrationResidenceBlockageSinkTauSimpleFoamV4FE
  *  go to applications/solvers/incompressible/wDissipationTempConcentrationResidenceBlockageSinkTauSimpleFoamV4FE
  *  run wmake
* also you'll need to compile the DoubleGermanCarreauConcentration material model (among others) from /src/transportModels/incompressible using wmake libso
  *  go to '/src/transportModels/incompressible/'
  *  run 'wmake libso'
* then head for the tutorialCoextrusion folder
* having arrived there, follow the local readMe instructions


