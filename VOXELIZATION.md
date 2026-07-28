# CitySim fixed-cell voxelization (V1)

This module converts CitySim's polygon scene into a compact, LEGO-like
Cartesian material grid and attaches a neutral log-law wind profile with V1
tree-canopy attenuation. It is a static preprocessing step; it does not yet
solve pollutant transport or building-induced flow.

## XML activation

Voxelization is opt-in. Add one `Voxelization` element directly under the
`CitySim` root, at the same level as `Simulation` and `District`:

```xml
<CitySim>
    <Simulation beginMonth="1" endMonth="1" beginDay="1" endDay="1"/>
    <Climate location="weather.cli"/>
    <Voxelization deltaX="2" deltaY="2" deltaZ="1"
                  roughness="0.3"/>
    <District>
        <Trees>
            <Tree id="1" leafAreaIndex="1" leafDistance="1"
                  attenuation="0.3">
                ...
            </Tree>
        </Trees>
    </District>
</CitySim>
```

All three dimensions and `roughness` are required. Dimensions and roughness
are expressed in metres and must be finite and strictly positive. Roughness
must also be below the fixed 10 m wind reference height. Each Tree may
independently override its optional dimensionless `attenuation`, whose default
is `0.5`. A readable climate file is required because its wind speed and
direction define the steady wind field. If `Voxelization` is absent, a normal
CitySim run creates no voxel file. Invalid or duplicate elements stop
processing with an error.

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

This mode requires the XML `Voxelization` and `Climate` elements and stops
after writing the grid. It is the fastest way to inspect the current
preprocessing result. The optional output argument changes only the file
location; it does not override the XML parameters.

The VTK file is a binary legacy `STRUCTURED_POINTS` data set with these cell
arrays:

- `material`: one unsigned byte per cell.
- `wind_speed`: speed magnitude in m/s as a big-endian float.
- `wind_velocity`: `(u,v,w)` in m/s as three big-endian floats.

| Value | Material | Pollutant unknown later? |
|---:|---|---:|
| 0 | `AIR` | yes |
| 1 | `BUILDING` | no |
| 2 | `GROUND` | no |
| 3 | `TREE` | yes |

Keeping building and ground distinct lets a later finite-volume solver obtain
the boundary condition directly from the material on the other side of each
air-cell face. Wind values are also written in solid cells to preserve regular
cell arrays; a later solver must ignore velocity unknowns in `BUILDING` and
`GROUND`.

## Public C++ API

`voxel.h` provides:

- `VoxelizationOptions`: resolution, repair, threading and allocation options.
- `SceneVoxelizer::voxelize(...)`: static scene conversion.
- `VoxelGrid`: contiguous material storage, interpolated column ground heights,
  compact per-level wind values, counts and VTK export.

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

### `roughness`

Required aerodynamic roughness length `z0` in metres. It must satisfy
`0 < z0 < 10`. No default is silently selected because terrain roughness is a
user/model choice.

The climate `FF` value at the simulation `beginDay`, hour 1 is interpreted as
the wind speed at `zref = 10 m`. Following the Beamer V1 model, the speed at a
voxel-centre height `z` is:

```text
Ulog(z) = U10 * ln((z + z0)/z0) / ln((10 + z0)/z0)
```

Height is measured from the lower domain plane and therefore the unattenuated
profile depends only on voxel level `k`. This deliberately ignores local
terrain-height changes and building flow disturbance in V1.

Climate `DD` is treated as a meteorological direction: degrees clockwise from
North describing where wind comes from. With CitySim `x=East`, `y=North`:

```text
u = -U sin(DD)
v = -U cos(DD)
w = 0
```

### Tree `attenuation`

Optional dimensionless canopy-class factor `a_t` on each `Tree`, with
`0 < a_t < 1`. Default: `0.5`. The Beamer gives `0.7`, `0.5` and `0.3` as
example sparse, medium and dense canopy classes respectively.

The cell multiplier and stored velocity are:

```text
aP = that Tree's attenuation  for its TREE cells
aP = 1                for every other cell
uP = aP * ulog(zP)
```

If multiple Tree crowns occupy one voxel, V1 uses the smallest attenuation
factor (the strongest obstruction). It does not multiply factors, which would
make the result depend on duplicate or overlapping Tree geometry.

Thus only TREE cells are slowed. Air immediately downstream returns to the
diagnostic log profile: V1 has no cumulative canopy wake, canopy drag,
pollutant uptake or vegetation deposition. No additional canopy drag or
leaf-density inputs are introduced at this stage. CitySim's existing integer
`leafAreaIndex` remains solely the legacy control that generates leaf layers;
the pollutant model does not reinterpret or estimate a physical quantity from
it.

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
- Wind uses only one speed and one three-component vector per z level.
  Per-tree factors are retained as sorted sparse `(cell index, attenuation)`
  arrays only for TREE cells. Scaling is applied while the expanded cell
  arrays are streamed row-by-row to VTK.
- Geometry is triangulated only once.
- Triangle tests are restricted to triangle bounding boxes.
- Ground queries use coarse spatial bins.
- Missing terrain reconstruction is linear in the number of horizontal cells.
- Recognized buildings use small local masks rather than one mask per building
  over the complete domain.
- Buildings are processed in OpenMP dynamic scheduling.
- Large column and material passes use OpenMP static scheduling.
- The material VTK payload is written contiguously; big-endian wind rows are
  generated and streamed without allocating full per-cell float arrays.
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

Following the Beamer, the later solver should use one shared harmonic
attenuation on each fluid-fluid face:

```text
af = 2*aP*aN/(aP+aN)
```

This modifies only advective flux. It is not a pollutant sink.

This avoids storing six face labels per cell and keeps preprocessing memory
small.

## Regression test

```text
make test-voxel
```

The test scene contains a building, two separated ground polygons, two leaf
levels and a thin trunk. It verifies binary VTK dimensions/materials, the
big-endian scalar/vector wind arrays, the shifted Beamer log law and exact
TREE-only canopy attenuation, including an explicit per-tree value and the
`0.5` default. It also checks ground gap interpolation, inclusive vertical
TREE filling, and that the trunk does not create solid or TREE cells below the
canopy. Missing/invalid voxel dimensions, roughness or climate input, and
invalid Tree attenuation create no VTK file and return an error.
Additional cases remove all explicit ground beneath a building and trees, and
exercise a tree-only scene; both must produce a continuous solid terrain layer
across the full domain.
