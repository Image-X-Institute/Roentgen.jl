```@setup abc
using Roentgen
```

# Types of Dose Positions

Dose positions are positions in the patient volume where dose is computed.
While this could span the whole patient volume, for some applications a smaller dose volume is preferred which target various regions or structures in the body.

A number of dose position types are provided.
There are generally two ways of creating dose positions: through the default constructor, or using a bounding object.


## Dose Grids

Dose grids are a subtype of dose positions where each position is located on a Cartesian grid.
Each axis need not be uniform, they can take any vector of positions. 
These are useful as they reduce memory usage, and allow for easy export to visualisation software.

### DoseGrid

`DoseGrid` is a simple Cartesian grid where dose is computed on every point in the box.

They are constructed by supplying grid axes:
```@repl abc
axes = (-10.:5.:10, -10.:5.:10, -10.:5.:10)
pos = DoseGrid(axes);
```

When using a bounding object, the dose grid dimensions are set to encompass the whole bounding object.
For example, a dose grid bounded by a cylinder will produce the following:
```@repl abc
pos = DoseGrid(5., CylinderBounds(10., 10.));
getaxes(pos)
```

It supports both linear and Cartesian indexing:
```@repl abc
pos[1]
pos[1, 2, 3]
```

### Masked Dose Grid

Masked dose grids are Cartesian grids where points can be masked such that they are not used when iterating through the grid.
This is an efficient implementation of sparse-like grids which maintain their 3D grid structure while avoiding computations on unnecessary positions.

Masked dose grids can be constructed by suppling grid axes and a vector of `CartesianIndex` where the points are enabled:
```@repl abc
indices = [CartesianIndex(1, 1, 1), CartesianIndex(2, 4, 1)]
pos = DoseGridMasked(axes, indices);
```
In this example, only positions at indices `[1, 1, 1]` and `[2, 4, 1]` are not masked and used in iteration:
```@repl abc
for i in eachindex(pos)
    @show i, pos[i]
end
```

When using a bounding object, the dose grid dimensions are set to encompass the whole bounding object, but points outside the bounding object are masked.
For example, a masked dose grid bounded by a cylinder will produce the following:
```@repl abc
pos = DoseGridMasked(5., CylinderBounds(10., 10.));
getaxes(pos)
length(pos)
```
While the axes still encompass the cylinder, the number of points (`length`) are less (15<5^3).

Linear indexing is supported, and Cartesian indexing will work if the index is not masked.
For example the following will work,
```@repl abc
pos[2]
pos[1,1,1]
```
But this will throw an `BoundsError` error,
```@repl abc
pos[1,2,1]
```
