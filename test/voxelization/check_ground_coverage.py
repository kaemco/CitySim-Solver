#!/usr/bin/env python3
"""Checks that an inferred terrain creates a continuous solid bottom layer."""

import pathlib
import sys

from check_voxels import index, read_vtk


def main():
    dimensions, values = read_vtk(pathlib.Path(sys.argv[1]))
    nx, ny, _ = dimensions
    bottom = [values[index(dimensions, i, j, 0)]
              for j in range(ny) for i in range(nx)]

    # A building can occupy a terrain-surface voxel, so both solid labels are
    # acceptable. AIR and TREE at the domain bottom indicate missing terrain.
    assert set(bottom) <= {1, 2}, set(bottom)
    assert 2 in bottom
    assert set(values) <= {0, 1, 2, 3}

    if len(sys.argv) > 2 and sys.argv[2] == "--tree-only":
        assert 1 not in values, "The ignored trunk became BUILDING"
        assert 0 in values and 2 in values and 3 in values

    print(
        "Complete terrain coverage passed:",
        dimensions,
        {"BUILDING_BOTTOM": bottom.count(1), "GROUND_BOTTOM": bottom.count(2)},
    )


if __name__ == "__main__":
    main()
