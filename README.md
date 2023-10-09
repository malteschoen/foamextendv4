# foamextendv4

# what is this?
we have solvers for extrusion operation (stationary, incompressible) taking into account both solids and fluids. More precisely, we have ...
1. ... immersed body solid mechanics taking into account the fluid forces (pressure, shear stress). This is a two-step, non-coupled approach - i.e. after a flow calculation there is no feedback of the geometry deformation into the fluid simulation. 
2. ... immersed body fluid flow in which fluid and solid domain share a heat/temperature equation. The funky stuff is the addition of heat sinks.
3. ... a coextrusion model based on concentration tracking - depending on concentration, different values of the Carreau-WLF-model are used in each cell

# i wanna try solid mechanics

# todo
organise nice tutorials for
1. 	just solid mechanics
2.  just heatsink
3.  just coextrusion
4.  a combination of the above if possible
