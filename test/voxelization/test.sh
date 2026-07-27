#!/bin/bash

set -e

solver="${1:-./CitySim}"
temporary_directory="$(mktemp -d)"
trap 'rm -rf "$temporary_directory"' EXIT

enabled_input="$temporary_directory/enabled.xml"
enabled_output="$temporary_directory/enabled.vtk"
cp test/voxelization/scene.xml "$enabled_input"
"$solver" --voxelize "$enabled_input" "$enabled_output"
python3 test/voxelization/check_voxels.py "$enabled_output"

# Remove every explicit Ground polygon. The non-simulated building floor must
# seed the terrain, which is then reconstructed below the building and trees.
no_ground_input="$temporary_directory/no_ground.xml"
no_ground_output="$temporary_directory/no_ground.vtk"
sed '/<GroundSurface /,/<\/GroundSurface>/d' \
    test/voxelization/scene.xml > "$no_ground_input"
"$solver" --voxelize "$no_ground_input" "$no_ground_output"
python3 test/voxelization/check_ground_coverage.py "$no_ground_output"

# A scene containing only leaves and a trunk must also obtain terrain. The trunk
# supplies only its base elevation and must remain absent from the material grid.
tree_output="$temporary_directory/tree_no_ground.vtk"
"$solver" --voxelize test/voxelization/tree_no_ground.xml "$tree_output"
python3 test/voxelization/check_ground_coverage.py \
    "$tree_output" --tree-only

missing_input="$temporary_directory/missing.xml"
missing_output="$temporary_directory/missing.vtk"
sed '/<Voxelization /d' test/voxelization/scene.xml > "$missing_input"
if "$solver" --voxelize "$missing_input" "$missing_output"; then
    echo "Expected a missing Voxelization element to be rejected." >&2
    exit 1
fi
if [ -e "$missing_output" ]; then
    echo "Voxel output was created without a Voxelization element." >&2
    exit 1
fi

invalid_input="$temporary_directory/invalid.xml"
invalid_output="$temporary_directory/invalid.vtk"
sed 's/deltaX="1"/deltaX="0"/' test/voxelization/scene.xml > "$invalid_input"
if "$solver" --voxelize "$invalid_input" "$invalid_output"; then
    echo "Expected a non-positive voxel dimension to be rejected." >&2
    exit 1
fi
if [ -e "$invalid_output" ]; then
    echo "Voxel output was created with an invalid voxel dimension." >&2
    exit 1
fi

# The ordinary simulation entry point must also honour the opt-in tag.
normal_directory="$temporary_directory/normal"
mkdir -p "$normal_directory" "$temporary_directory/simulations/ClimateFile"
cp test/voxelization/scene.xml "$normal_directory/enabled.xml"
cp test/simulations/ClimateFile/DRYCOLD3.cli \
   "$temporary_directory/simulations/ClimateFile/DRYCOLD3.cli"
"$solver" "$normal_directory/enabled.xml"
python3 test/voxelization/check_voxels.py "$normal_directory/enabled_voxels.vtk"

# Without the tag, an ordinary simulation must not create a voxel file.
sed '/<Voxelization /d' test/voxelization/scene.xml > "$normal_directory/missing.xml"
"$solver" "$normal_directory/missing.xml"
if [ -e "$normal_directory/missing_voxels.vtk" ]; then
    echo "Ordinary simulation created voxel output without a Voxelization element." >&2
    exit 1
fi
