```@setup abc
using Roentgen
```

# Generating Fluence Maps

In order to compute dose, an idealised fluence distribution needs to be created.

This is usually done by generating fluence maps on a set of bixels from beam-limiting devices, such as the jaws or MLC.

## Bixels

Bixels are 2D elements denoting a rectangular section of the isoplane.
Their position and width is in the IEC BLD coordinate system, scaled to the isoplane.

Bixels are constructed

## Beam-Limiting Devices

Beam-limiting devices contain information on the position and shape of the aperture.
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
The leaves move in the x direction.

They are constructed by specifying leaf positions in a `2xn` matrix and leaf edges in an `n+1` length vector,
```@repl abc
mlcx = [-10.  -8. -9.
         10.   5.  0.];
mlcy = [-5., 0., 5., 10.];
mlc = MultiLeafCollimator(mlcx, mlcy);
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



## Fluence

