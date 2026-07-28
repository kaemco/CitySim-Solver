#!/bin/bash

set -e
export PYTHONDONTWRITEBYTECODE=1

solver="${1:-./CitySim}"
temporary_directory="$(mktemp -d)"
trap 'rm -rf "$temporary_directory"' EXIT

mkdir -p "$temporary_directory/simulations/ClimateFile"
cp test/simulations/ClimateFile/DRYCOLD3.cli \
   "$temporary_directory/simulations/ClimateFile/DRYCOLD3.cli"

voxel_directory="$temporary_directory/voxel"
mkdir -p "$voxel_directory"

enabled_input="$voxel_directory/enabled.xml"
enabled_output="$voxel_directory/enabled.vtk"
cp test/voxelization/scene.xml "$enabled_input"
"$solver" --voxelize "$enabled_input" "$enabled_output"
python3 test/voxelization/check_voxels.py "$enabled_output" 0.3

# Remove every explicit Ground polygon. The non-simulated building floor must
# seed the terrain, which is then reconstructed below the building and trees.
no_ground_input="$voxel_directory/no_ground.xml"
no_ground_output="$voxel_directory/no_ground.vtk"
sed '/<GroundSurface /,/<\/GroundSurface>/d' \
    test/voxelization/scene.xml > "$no_ground_input"
"$solver" --voxelize "$no_ground_input" "$no_ground_output"
python3 test/voxelization/check_ground_coverage.py "$no_ground_output"

# A scene containing only leaves and a trunk must also obtain terrain. The trunk
# supplies only its base elevation and must remain absent from the material grid.
tree_input="$voxel_directory/tree_no_ground.xml"
tree_output="$voxel_directory/tree_no_ground.vtk"
cp test/voxelization/tree_no_ground.xml "$tree_input"
"$solver" --voxelize "$tree_input" "$tree_output"
python3 test/voxelization/check_ground_coverage.py \
    "$tree_output" --tree-only

# Different Tree elements in one scene must preserve their own factors.
two_tree_input="$voxel_directory/two_tree_attenuations.xml"
two_tree_output="$voxel_directory/two_tree_attenuations.vtk"
cp test/voxelization/two_tree_attenuations.xml "$two_tree_input"
"$solver" --voxelize "$two_tree_input" "$two_tree_output"
python3 test/voxelization/check_two_tree_attenuations.py "$two_tree_output"

missing_input="$voxel_directory/missing.xml"
missing_output="$voxel_directory/missing.vtk"
sed '/<Voxelization /d' test/voxelization/scene.xml > "$missing_input"
if "$solver" --voxelize "$missing_input" "$missing_output"; then
    echo "Expected a missing Voxelization element to be rejected." >&2
    exit 1
fi
if [ -e "$missing_output" ]; then
    echo "Voxel output was created without a Voxelization element." >&2
    exit 1
fi

invalid_input="$voxel_directory/invalid.xml"
invalid_output="$voxel_directory/invalid.vtk"
sed 's/deltaX="1"/deltaX="0"/' test/voxelization/scene.xml > "$invalid_input"
if "$solver" --voxelize "$invalid_input" "$invalid_output"; then
    echo "Expected a non-positive voxel dimension to be rejected." >&2
    exit 1
fi
if [ -e "$invalid_output" ]; then
    echo "Voxel output was created with an invalid voxel dimension." >&2
    exit 1
fi

missing_roughness_input="$voxel_directory/missing_roughness.xml"
missing_roughness_output="$voxel_directory/missing_roughness.vtk"
sed 's/ roughness="0.1"//' \
    test/voxelization/scene.xml > "$missing_roughness_input"
if "$solver" --voxelize "$missing_roughness_input" "$missing_roughness_output"; then
    echo "Expected a missing roughness length to be rejected." >&2
    exit 1
fi
if [ -e "$missing_roughness_output" ]; then
    echo "Voxel output was created without a roughness length." >&2
    exit 1
fi

invalid_roughness_input="$voxel_directory/invalid_roughness.xml"
invalid_roughness_output="$voxel_directory/invalid_roughness.vtk"
sed 's/roughness="0.1"/roughness="10"/' \
    test/voxelization/scene.xml > "$invalid_roughness_input"
if "$solver" --voxelize "$invalid_roughness_input" "$invalid_roughness_output"; then
    echo "Expected roughness at the reference height to be rejected." >&2
    exit 1
fi
if [ -e "$invalid_roughness_output" ]; then
    echo "Voxel output was created with an invalid roughness length." >&2
    exit 1
fi

default_tree_attenuation_input="$voxel_directory/default_tree_attenuation.xml"
default_tree_attenuation_output="$voxel_directory/default_tree_attenuation.vtk"
sed 's/ attenuation="0.3"//' \
    test/voxelization/scene.xml > "$default_tree_attenuation_input"
"$solver" --voxelize \
    "$default_tree_attenuation_input" "$default_tree_attenuation_output"
python3 test/voxelization/check_voxels.py "$default_tree_attenuation_output"

zero_tree_attenuation_input="$voxel_directory/zero_tree_attenuation.xml"
zero_tree_attenuation_output="$voxel_directory/zero_tree_attenuation.vtk"
sed 's/attenuation="0.3"/attenuation="0"/' \
    test/voxelization/scene.xml > "$zero_tree_attenuation_input"
if "$solver" --voxelize \
        "$zero_tree_attenuation_input" "$zero_tree_attenuation_output"; then
    echo "Expected zero tree attenuation to be rejected." >&2
    exit 1
fi
if [ -e "$zero_tree_attenuation_output" ]; then
    echo "Voxel output was created with zero tree attenuation." >&2
    exit 1
fi

unit_tree_attenuation_input="$voxel_directory/unit_tree_attenuation.xml"
unit_tree_attenuation_output="$voxel_directory/unit_tree_attenuation.vtk"
sed 's/attenuation="0.3"/attenuation="1"/' \
    test/voxelization/scene.xml > "$unit_tree_attenuation_input"
if "$solver" --voxelize \
        "$unit_tree_attenuation_input" "$unit_tree_attenuation_output"; then
    echo "Expected unit tree attenuation to be rejected." >&2
    exit 1
fi
if [ -e "$unit_tree_attenuation_output" ]; then
    echo "Voxel output was created with unit tree attenuation." >&2
    exit 1
fi

missing_climate_input="$voxel_directory/missing_climate.xml"
missing_climate_output="$voxel_directory/missing_climate.vtk"
sed '/<Climate /d' test/voxelization/scene.xml > "$missing_climate_input"
if "$solver" --voxelize "$missing_climate_input" "$missing_climate_output"; then
    echo "Expected a missing Climate element to be rejected." >&2
    exit 1
fi
if [ -e "$missing_climate_output" ]; then
    echo "Voxel output was created without climate data." >&2
    exit 1
fi

# The ordinary simulation entry point must also honour the opt-in tag.
normal_directory="$temporary_directory/normal"
mkdir -p "$normal_directory"
cp test/voxelization/scene.xml "$normal_directory/enabled.xml"
"$solver" "$normal_directory/enabled.xml"
python3 test/voxelization/check_voxels.py \
    "$normal_directory/enabled_voxels.vtk" 0.3

# Without the tag, an ordinary simulation must not create a voxel file.
sed '/<Voxelization /d' test/voxelization/scene.xml > "$normal_directory/missing.xml"
"$solver" "$normal_directory/missing.xml"
if [ -e "$normal_directory/missing_voxels.vtk" ]; then
    echo "Ordinary simulation created voxel output without a Voxelization element." >&2
    exit 1
fi
