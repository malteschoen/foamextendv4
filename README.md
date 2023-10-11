# foamextendv4
this should run with foam-extend version 4
todo: check performance in fe-v4.1 and fe-v5.x

# what is this?
we have solvers for extrusion operation (stationary, incompressible) taking into account both solids and fluids. More precisely, we have ...
1. ... immersed body solid mechanics taking into account the fluid forces (pressure, shear stress). This is a two-step, non-coupled approach - i.e. after a flow calculation there is no feedback of the geometry deformation into the fluid simulation. 
2. ... immersed body fluid flow in which fluid and solid domain share a heat/temperature equation. The funky stuff is the addition of heat sinks.
3. ... a coextrusion model based on concentration tracking - depending on concentration, different values of the Carreau-WLF-model are used in each cell

# i just wanna try solid mechanics
* splendid!
* before you go ahead, you'll need to compile wDissipationTempConcentrationResidenceBlockageSinkTauSimpleFoamV4FE from applications/solvers/incompressible using wmake.
* also, you will need to compile /applications/solvers/solidMechanics/immersedTauFastenerElasticSolidFoam using wmake
* having arrived there, follow the local readMe instructions


# i just wanna try the heat-sink stuff
* excellent!
* before you go ahead, you'll need to compile wDissipationTempConcentrationResidenceBlockageSinkTauSimpleFoamV4FE from applications/solvers/incompressible using wmake.
* then head for the tutorialHeatSink folder
* having arrived there, follow the local readMe instructions

# i just wanna try coextrusion
* magnificient!
* before you go ahead, you'll need to compile wDissipationTempConcentrationResidenceBlockageSinkTauSimpleFoamV4FE from applications/solvers/incompressible using wmake.
* also you'll need to compile the DoubleGermanCarreau material model from \src\transportModels\incompressible using wmake libso
* then head for the tutorialCoextrusion folder
* having arrived there, follow the local readMe instructions

# todo
organise nice tutorials for
1. 	just solid mechanics
2.  just heatsink
3.  just coextrusion
4.  a combination of the above if possible
