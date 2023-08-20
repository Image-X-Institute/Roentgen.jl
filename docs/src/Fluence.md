```@setup abc
using Roentgen
```

# Generating Fluence Maps

In order to compute dose, an idealised fluence distribution needs to be created.

This is usually done by generating fluence maps on a set of bixels from beam-limiting devices, such as the jaws or MLC.

## Beam-Limiting Devices

Beam-limiting devices contain information on the position and shape of the aperture.
These typically relate to physical devices on the treatment machine, such as a Multi-Leaf Collimator.
Each have specialised methods by which they interact with bixels and assign an amount of fluence going through the bixel.

### Jaws

Jaws are rectangular field shapes, which are simply constructed by specifying the position of the left, right, top and bottom jaws, or a fieldsize:
```@repl abc
jaws = Jaws(-10., 10., -15., 15.)
jaws = Jaws(100.)
```

They have methods for accessing the positions, and the area,
```@repl abc
getx(jaws)
gety(jaws)
getarea(jaws)
```

### Multi-Leaf Collimator

Multi-Leac Collimator (MLC) contain the position of the leaves (in the x direction), and the position of the leaf edges (the y direction).

They are constructed by specifying leaf positions in a `2xn` matrix and leaf edges in an `n+1` length vector,
```@repl abc
mlcx = [-10. -12. -23.
         10.   4. -10.];
mlcy = -5.:5.:10;
mlc = MultiLeafCollimator(mlcx, mlcy)
```
They can also be constructed by just specifying the edges, which sets the leaf positions to zero,
```julia
MultiLeafCollimator(mlcy)
```

Indexing and iterating through the MLC will produce a subset of the MLC, either returning `Jaws` or another `MultiLeafCollimator`,
```@repl abc
mlc[1]
mlc[2:3]
@view mlc[1:2]
```

Leaf positions can be accessed using the `getpositions`,
```@repl abc
getpositions(mlc)
getpositions(mlc, 1)
```

Edge positions can be accessed using `getedges`,
```@repl abc
getedges(mlc)
getedges(mlc, 1)
getedges(mlc, 1:2)
```

Leaf positions can be set with `setpositions!` and closed with `closeleaves!`
```@repl abc
closeleaves!(mlc)
mlc[1] = -5., 25
mlc
mlcx = [-10. -12. -23.
         10.   4. -10.];
setpositions!(mlc, mlcx)
```

## Bixels

Bixels are 2D elements denoting a rectangular section of the isoplane.
Their position and width is in the IEC BLD coordinate system, scaled to the isoplane.

Bixels are constructed by specifing their position and width,
```@repl abc
bixel = Bixel(0., 0., 1., 1.)
```

A collection of bixels is used to build a fluence map.
Most functions will take a list of bixels (`AbstractVector{Bixel}`).
These can either be constructed using the `Bixel` constructor, or using one of the following bixel methods.

### Bixel Grids

Bixel grids store bixels in a rectilinear grid, and behave as `Matrix{Bixel}`.
They can be constructed by supplying a vector or range of bixel edge positions,
```@repl abc
bixels = BixelGrid(-10.:5.:10, -10.:5.:10.)
```
Bixel grids can also be constructed using the information from beam-limiting devices,
```@repl abc
jaws = Jaws([-5., 5.], [-10., 10.])
bixels = BixelGrid(jaws, 5., 5.)
```

They behave like matrices: they can be indexed (both linear and Cartesian indexing), and iterate through each bixel,
```@repl abc
bixels[3]
bixels[1, 2]
```

They also have specialised methods for accessing grid axes,
```@repl abc
getaxes(bixels)
getaxes(bixels, 1)
```

## Fluence

