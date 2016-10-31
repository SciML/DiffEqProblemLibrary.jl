
###Example Meshes
pkgdir = dirname(@__FILE__)
meshes_location = "premade_meshes.jld"
"`meshExample_bunny()` : Returns a 3D SimpleMesh."
meshExample_bunny() = load("$pkgdir/$meshes_location","bunny")

"`meshExample_flowpastcylindermesh()` : Returns a 2D SimpleMesh."
meshExample_flowpastcylindermesh() = load("$pkgdir/$meshes_location","flowpastcylindermesh")

"`meshExample_lakemesh()` : Returns a 2D SimpleMesh."
meshExample_lakemesh() = load("$pkgdir/$meshes_location","lakemesh")

"`meshExample_Lshapemesh()` : Returns a 2D SimpleMesh."
meshExample_Lshapemesh() = load("$pkgdir/$meshes_location","Lshapemesh")

"`meshExample_Lshapeunstructure()` : Returns a 2D SimpleMesh."
meshExample_Lshapeunstructure() = load("$pkgdir/$meshes_location","Lshapeunstructure")

"`meshExample_oilpump()` : Returns a 3D SimpleMesh."
meshExample_oilpump() = load("$pkgdir/$meshes_location","oilpump")

"`meshExample_wavymesh()` : Returns a 2D SimpleMesh."
meshExample_wavymesh() = load("$pkgdir/$meshes_location","wavymesh")

"`meshExample_wavyperturbmesh()` : Returns a 3D SimpleMesh."
meshExample_wavyperturbmesh() = load("$pkgdir/$meshes_location","wavyperturbmesh")
