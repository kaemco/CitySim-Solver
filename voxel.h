#ifndef VOXEL_H
#define VOXEL_H

#include <array>
#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <limits>
#include <string>
#include <vector>

class XmlScene;

/**
 * Material values intentionally fit in one byte.  BUILDING and GROUND remain
 * distinct so pollutant boundary faces can later apply different conditions.
 */
enum class VoxelMaterial : std::uint8_t {
    AIR = 0,
    BUILDING = 1,
    GROUND = 2,
    TREE = 3
};

struct VoxelizationOptions {
    double deltaX = 2.;
    double deltaY = 2.;
    double deltaZ = 2.;

    /**
     * Number of one-cell morphological closing passes applied independently
     * to every Building.  This repairs grid-scale cracks without modifying the
     * original CitySim surfaces.  Set to zero for maximum speed/no repair.
     */
    unsigned int buildingGapClosingPasses = 1;

    /**
     * Coarse horizontal bin size used to accelerate ground-triangle queries.
     */
    unsigned int groundTileSize = 16;

    /**
     * 0 uses the OpenMP runtime default.  Values greater than zero request an
     * explicit number of threads for parallel loops.
     */
    int threads = 0;

    /**
     * Generic District surfaces include shading surfaces and the surfaces of
     * Building elements declared Simulate="false".  In V1 they are treated as
     * thin BUILDING obstacles so the latter are not lost.
     */
    bool includeDistrictSurfacesAsSolids = true;

    /**
     * Allocation guard.  The compact material grid itself uses one byte/cell,
     * but geometry processing needs temporary storage as well.
     */
    std::size_t maximumCellCount = 500000000u;
};

struct VoxelCounts {
    std::size_t air = 0;
    std::size_t building = 0;
    std::size_t ground = 0;
    std::size_t tree = 0;
};

class VoxelGrid {
public:
    VoxelGrid();
    VoxelGrid(double originX, double originY, double originZ,
              double deltaX, double deltaY, double deltaZ,
              std::size_t countX, std::size_t countY, std::size_t countZ);

    std::size_t countX() const { return countX_; }
    std::size_t countY() const { return countY_; }
    std::size_t countZ() const { return countZ_; }
    std::size_t cellCount() const { return materials_.size(); }

    double originX() const { return originX_; }
    double originY() const { return originY_; }
    double originZ() const { return originZ_; }
    double deltaX() const { return deltaX_; }
    double deltaY() const { return deltaY_; }
    double deltaZ() const { return deltaZ_; }
    double maximumX() const { return originX_ + deltaX_*countX_; }
    double maximumY() const { return originY_ + deltaY_*countY_; }
    double maximumZ() const { return originZ_ + deltaZ_*countZ_; }

    std::size_t index(std::size_t i, std::size_t j, std::size_t k) const {
        return (k*countY_ + j)*countX_ + i;
    }

    VoxelMaterial material(std::size_t i, std::size_t j, std::size_t k) const {
        return static_cast<VoxelMaterial>(materials_.at(index(i,j,k)));
    }

    const std::vector<std::uint8_t>& materials() const { return materials_; }
    std::vector<std::uint8_t>& materials() { return materials_; }

    /**
     * Reconstructed ground elevation for each (i,j) column. A completed
     * SceneVoxelizer result has a finite value in every horizontal column.
     */
    float groundHeight(std::size_t i, std::size_t j) const {
        return groundHeights_.at(j*countX_ + i);
    }
    const std::vector<float>& groundHeights() const { return groundHeights_; }
    std::vector<float>& groundHeights() { return groundHeights_; }

    VoxelCounts counts() const;

    /**
     * Fast legacy VTK STRUCTURED_POINTS writer.  Materials are unsigned bytes:
     * 0=AIR, 1=BUILDING, 2=GROUND, 3=TREE.
     */
    void writeLegacyVtk(const std::string& fileName) const;

private:
    double originX_ = 0.;
    double originY_ = 0.;
    double originZ_ = 0.;
    double deltaX_ = 1.;
    double deltaY_ = 1.;
    double deltaZ_ = 1.;
    std::size_t countX_ = 0;
    std::size_t countY_ = 0;
    std::size_t countZ_ = 0;
    std::vector<std::uint8_t> materials_;
    std::vector<float> groundHeights_;
};

class SceneVoxelizer {
public:
    /**
     * Converts the physical CitySim geometry to a fixed-cell, LEGO-like grid.
     * The grid covers the geometry AABB without physical-model padding; only
     * the upper grid face may extend by less than one cell to keep fixed sizes.
     */
    static VoxelGrid voxelize(XmlScene& scene,
                              const VoxelizationOptions& options,
                              std::ostream* logStream = nullptr);
};

#endif
