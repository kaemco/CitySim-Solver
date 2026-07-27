# CitySim fixed-cell voxelization (V1)

This module converts CitySim's polygon scene into a compact, LEGO-like
Cartesian material grid. It is a static geometry preprocessing step; it does
not yet solve wind or pollutant transport.

## XML activation

Voxelization is opt-in. Add one `Voxelization` element directly under the
`CitySim` root, at the same level as `Simulation` and `District`:

```xml
<CitySim>
    <Simulation beginMonth="1" endMonth="1" beginDay="1" endDay="1"/>
    <Voxelization deltaX="2" deltaY="2" deltaZ="1"/>
    <District>
        ...
    </District>
</CitySim>
```

All three dimensions are required, expressed in metres, and must be finite and
strictly positive. If the element is absent, a normal CitySim run performs no
voxelization and creates no voxel file. Invalid or duplicate elements stop
processing with an error instead of silently using a default.

During a normal run:

```text
CitySim district.xml
```

CitySim writes `<scene>_voxels.vtk` before starting the existing simulations,
then continues normally.

## Geometry-only verification command

```text
CitySim --voxelize scene.xml [output.vtk]
```

This mode also requires the XML `Voxelization` element, but does not require a
climate file and stops after writing the grid. It is the fastest way to inspect
the current preprocessing result. The optional output argument changes only the
file location; it does not override the XML cell dimensions.

The VTK file is a binary legacy `STRUCTURED_POINTS` data set with one unsigned
byte named `material` per cell:

| Value | Material | Pollutant unknown later? |
|---:|---|---:|
| 0 | `AIR` | yes |
| 1 | `BUILDING` | no |
| 2 | `GROUND` | no |
| 3 | `TREE` | yes |

Keeping building and ground distinct lets a later finite-volume solver obtain
the boundary condition directly from the material on the other side of each
air-cell face.

## Public C++ API

`voxel.h` provides:

- `VoxelizationOptions`: resolution, repair, threading and allocation options.
- `SceneVoxelizer::voxelize(...)`: static scene conversion.
- `VoxelGrid`: contiguous material storage, interpolated column ground heights,
  counts and VTK export.

Materials are stored with `x` as the fastest-changing index:

```text
index = (k * countY + j) * countX + i
```

This layout gives contiguous horizontal rows and is suitable for later
six-neighbour finite-volume stencils.

## Parameters

### `deltaX`, `deltaY`, `deltaZ`

Required positive, finite voxel dimensions in metres, read from the
`CitySim/Voxelization` element. Cells are always the requested size; there are
no smaller remainder cells.

### `threads`

The XML interface currently uses `0`, meaning the OpenMP runtime default. A
positive C++ API value is used by the parallel
ground queries, per-building voxelization, large material loops and tree-column
filling. A build without OpenMP remains functional and executes serially.

### `buildingGapClosingPasses`

Default: `1`.

Each pass performs a six-neighbour dilation followed by erosion on each
CitySim `Building` independently. Original shell cells are always retained.
This repairs cracks at approximately one-voxel scale without joining different
recognized buildings. Set it to `0` for maximum speed and no explicit repair.

Generic District surfaces cannot be grouped reliably by building and therefore
do not receive this repair.

### `groundTileSize`

C++ API default: `16` horizontal cells. Ground triangles are assigned to
coarse horizontal tiles so each grid column tests only nearby triangles.
Smaller tiles reduce triangle candidates but use more bin metadata.

### `includeDistrictSurfacesAsSolids`

C++ API default: `true`.

CitySim puts both shading polygons and the surfaces of
`Building Simulate="false"` into the District's generic surface collection.
V1 conservatively treats all of them as thin `BUILDING` obstacles. Setting
this option to `false` ignores both groups.

### `maximumCellCount`

C++ API default: `500,000,000`. This prevents accidental allocations larger
than the intended high-speed fixed grid. V1 also restricts the grid to fewer
than `2^32` cells because exterior flood-fill indices use 32 bits.

## Domain

The domain is the axis-aligned bounding box of the geometry actually used:

- Walls, roofs, floors and generic zone surfaces of simulated buildings.
- Generic District surfaces when enabled, including shading geometry and
  `Building Simulate="false"`. They therefore contribute to the scene bounding
  box as well as to the solid mask.
- Ground surfaces.
- Explicit and generated leaf surfaces.

Tree `Trunc` surfaces create no `TREE` or `BUILDING` cells and do not enlarge
the horizontal bounds. Their lowest vertex can lower the domain minimum in `z`
because it provides terrain-height evidence for a tree without explicit ground.

For each direction:

```text
count = max(1, ceil((sceneMaximum - sceneMinimum) / delta))
```

The grid origin is exactly the scene minimum. Its upper face can exceed the
scene maximum by less than one voxel so all cells remain fixed-size. There is
no physical-model padding.

Consequently, if the scene bounding box is exactly a closed building, the
result can contain only building cells and no outside air. Exterior air space
will require padding in the future transport domain.

## Geometry conversion

CitySim polygons are first triangulated with the existing
`RENIndexedFaceSetBuilder`. Degenerate triangles are discarded. Intersection
and grid-coordinate calculations use `double`, even though CitySim stores its
source vertices as `float`.

### Buildings and thin obstacles

Every triangle is tested against only the cells in its index-space bounding
box using a 13-axis triangle/AABB separating-axis test.

An intersected voxel becomes entirely `BUILDING`. The original continuous
surface is therefore replaced by a staircase boundary on the faces between
building and fluid voxels. Grid-aligned zero-thickness surfaces are assigned
deterministically to one candidate layer before the exact intersection test.

Recognized buildings are rasterized into separate local masks, optionally
closed, and merged into the global material grid. Generic District surfaces
are rasterized directly without morphological repair.

After all shells are present, a six-neighbour exterior flood fill starts from
all air cells on the grid boundary. Air not reachable from the boundary is
considered an enclosed building interior and becomes `BUILDING`.

This handles:

- Reversed or inconsistent polygon orientation.
- Small non-manifold cracks sealed by conservative voxel overlap or closing.
- Missing floors when another solid surface, normally ground, closes the shell.

It cannot infer a building interior through a large missing façade. Such a
shell remains connected to outside air.

### Ground

Ground is treated as a complete, single-valued terrain height field rather than
a closed input mesh. Every horizontal grid column receives a finite elevation,
even if the CitySim scene has no `Ground` polygons.

Terrain observations have this precedence:

1. An explicit `Ground` intersection is authoritative. For overlapping ground
   polygons, the uppermost elevation is retained.
2. In a column not covered by explicit ground, the lowest building `Floor`
   intersection is used. Floors from both simulated buildings and
   `Building Simulate="false"` are retained as support footprints.
3. A tree-trunk base is used as a sparse elevation observation when neither an
   explicit ground nor a building floor covers its column. The trunk remains
   absent from the material grid.
4. If the complete scene contains none of these observations, the scene
   minimum `z` is used as a deterministic flat-terrain fallback.

Triangle queries use the same coarse horizontal bins as the surface
voxelization. Missing columns between terrain observations are filled in two
linear-time passes:

1. Linear interpolation between known values along `x`, with nearest-value
   extrapolation at row ends.
2. Linear interpolation between valid rows along `y`, with nearest-row
   extrapolation at the domain edges.

This produces one finite terrain elevation for every column. All cells from the
domain bottom through the cell containing that elevation initially become
`GROUND`. A building shell can subsequently occupy the coincident terrain
surface cell; this represents a building resting directly on the ground rather
than an air gap.

The ground algorithm assumes terrain is single-valued in `z`. For overlapping
terrain observations, the precedence above is used. Caves and overhangs are not
represented.

The support-footprint idea follows CitySim's binary STL export for 3D printing,
which adds caps from building floors and explicit grounds. Unlike that exporter,
voxelization keeps their physical elevations and reconstructs the surrounding
terrain instead of projecting them to an artificial printing base.

### Trees

For every CitySim `Tree`, V1 uses:

- Explicit `Leaf` surfaces.
- CitySim-generated `subLeaves` produced from leaf-area-index layers.

`Trunc` surfaces are deliberately ignored: they create neither TREE nor solid
cells.

For each horizontal column touched by at least one leaf triangle, V1 records
the minimum and maximum leaf elevation for that same tree. Every air cell from
the voxel containing the minimum leaf through the voxel containing the maximum
leaf, inclusive, becomes `TREE`.

The horizontal tree footprint is therefore the union of the projected leaf
surfaces. Separate trees are processed independently, so vertically separated
overlapping trees do not create an artificial TREE column between them.

Solid precedence is:

```text
BUILDING or GROUND > TREE > AIR
```

Trees remain fluid cells for the later pollutant solve.

## Performance characteristics

- Geometry points use a fixed POD of three adjacent `double` values (24 bytes,
  no allocation). Triangles are held contiguously in `std::vector<Triangle>`.
  A `std::array<double,3>` would have the same layout, while a
  `std::vector<double>` per point would add a separate allocation and pointer,
  size and capacity metadata to every point.
- One byte per material cell in contiguous storage.
- Ground heights require one `float` per horizontal column.
- Geometry is triangulated only once.
- Triangle tests are restricted to triangle bounding boxes.
- Ground queries use coarse spatial bins.
- Missing terrain reconstruction is linear in the number of horizontal cells.
- Recognized buildings use small local masks rather than one mask per building
  over the complete domain.
- Buildings are processed in OpenMP dynamic scheduling.
- Large column and material passes use OpenMP static scheduling.
- The VTK payload is binary and written in one contiguous operation.
- Exterior connectivity is currently a serial depth-first traversal using
  temporary 32-bit indices; it is the main remaining non-parallel preprocessing
  stage.

## Boundary interfaces for the future solver

The voxelizer does not preallocate boundary-face objects. They are cheap to
derive while assembling the finite-volume stencil: inspect the six neighbours
of each `AIR` or `TREE` cell.

- Neighbour `BUILDING`: zero normal advective flux and zero diffusive flux.
- Neighbour `GROUND`: zero normal advective flux and a ground-deposition
  diffusive/Robin contribution.
- Neighbour `AIR` or `TREE`: ordinary fluid-fluid face.

This avoids storing six face labels per cell and keeps preprocessing memory
small.

## Regression test

```text
make test-voxel
```

The test scene contains a building, two separated ground polygons, two leaf
levels and a thin trunk. It verifies binary VTK dimensions/materials, ground
gap interpolation, inclusive vertical TREE filling, and that the trunk does
not create solid or TREE cells below the canopy. It also verifies that a
missing `Voxelization` element and a non-positive cell dimension create no VTK
file and return an error. Additional cases remove all explicit ground beneath a
building and trees, and exercise a tree-only scene; both must produce a
continuous solid terrain layer across the full domain.
