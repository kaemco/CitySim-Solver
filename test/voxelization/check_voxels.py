#!/usr/bin/env python3
"""Checks the binary legacy-VTK output of the voxelization regression scene."""

import pathlib
import sys


def read_vtk(path):
    with open(path, "rb") as stream:
        dimensions = None
        cell_count = None
        while True:
            line = stream.readline()
            if not line:
                raise RuntimeError("Unexpected end of VTK header")
            text = line.decode("ascii").strip()
            if text.startswith("DIMENSIONS "):
                points = tuple(map(int, text.split()[1:]))
                dimensions = tuple(value - 1 for value in points)
            elif text.startswith("CELL_DATA "):
                cell_count = int(text.split()[1])
            elif text == "LOOKUP_TABLE default":
                break
        if dimensions is None or cell_count is None:
            raise RuntimeError("Incomplete VTK header")
        values = stream.read(cell_count)
        if len(values) != cell_count:
            raise RuntimeError("Incomplete VTK material payload")
        return dimensions, values


def index(dimensions, i, j, k):
    nx, ny, _ = dimensions
    return (k * ny + j) * nx + i


def main():
    dimensions, values = read_vtk(pathlib.Path(sys.argv[1]))
    assert dimensions == (6, 4, 4), dimensions
    assert set(values) <= {0, 1, 2, 3}
    assert all(values.count(material) > 0 for material in (0, 1, 2, 3))

    # x=3.5 lies in the deliberate gap between the two ground polygons.
    assert values[index(dimensions, 3, 3, 0)] == 2

    # The trunk is ignored: below the lowest leaf this canopy column is air.
    assert values[index(dimensions, 5, 1, 1)] == 0

    # TREE cells form vertically continuous runs from the minimum leaf voxel
    # through the maximum leaf voxel, inclusive.
    nx, ny, nz = dimensions
    tree_columns = 0
    for j in range(ny):
        for i in range(nx):
            levels = [k for k in range(nz) if values[index(dimensions, i, j, k)] == 3]
            if levels:
                tree_columns += 1
                assert levels == list(range(min(levels), max(levels) + 1))
                assert min(levels) >= 2
    assert tree_columns > 0

    print(
        "Voxel regression passed:",
        dimensions,
        {name: values.count(value) for name, value in
         (("AIR", 0), ("BUILDING", 1), ("GROUND", 2), ("TREE", 3))},
    )


if __name__ == "__main__":
    main()
