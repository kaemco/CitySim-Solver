#ifndef VOXEL_H
#define VOXEL_H

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <limits>
#include <string>
#include <utility>
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
     * Aerodynamic roughness length z0 in metres for the neutral log-law wind
     * profile.  XML input requires this value explicitly.  It must be finite,
     * positive and smaller than the fixed 10 m climate reference height.
     */
    double roughnessLength = std::numeric_limits<double>::quiet_NaN();

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

    /**
     * The unattenuated neutral log-law needs only countZ values in memory.
     * Actual cell wind additionally applies its owning Tree's attenuation.
     */
    bool hasWindProfile() const {
        return windSpeedByLevel_.size()==countZ_
                && windVelocityByLevel_.size()==countZ_;
    }
    float logLawWindSpeed(std::size_t k) const {
        return windSpeedByLevel_.at(k);
    }
    const std::array<float,3>& logLawWindVelocity(std::size_t k) const {
        return windVelocityByLevel_.at(k);
    }
    float windAttenuation(std::size_t i, std::size_t j, std::size_t k) const {
        const std::uint32_t cell=static_cast<std::uint32_t>(index(i,j,k));
        const std::vector<std::uint32_t>::const_iterator found=
                std::lower_bound(treeWindCellIndices_.begin(),
                                 treeWindCellIndices_.end(),cell);
        if (found==treeWindCellIndices_.end() || *found!=cell) return 1.f;
        return treeWindAttenuations_[
                static_cast<std::size_t>(found-treeWindCellIndices_.begin())];
    }
    float windSpeed(std::size_t i, std::size_t j, std::size_t k) const {
        return logLawWindSpeed(k)*windAttenuation(i,j,k);
    }
    std::array<float,3> windVelocity(std::size_t i, std::size_t j,
                                     std::size_t k) const {
        const float attenuation=windAttenuation(i,j,k);
        const std::array<float,3>& base=logLawWindVelocity(k);
        return std::array<float,3>{{
            attenuation*base[0],attenuation*base[1],attenuation*base[2]
        }};
    }
    double referenceWindSpeed() const { return referenceWindSpeed_; }
    double referenceWindDirection() const { return referenceWindDirection_; }
    double roughnessLength() const { return roughnessLength_; }
    double referenceWindHeight() const { return referenceWindHeight_; }
    std::size_t attenuatedTreeCellCount() const {
        return treeWindCellIndices_.size();
    }

    void setLogLawWindProfile(double referenceSpeed,
                              double meteorologicalDirectionDegrees,
                              double roughnessLength,
                              double referenceHeight=10.);
    void setTreeWindAttenuations(
            std::vector<std::pair<std::uint32_t,float> > assignments);

    VoxelCounts counts() const;

    /**
     * Fast legacy VTK STRUCTURED_POINTS writer.  Materials are unsigned bytes
     * (0=AIR, 1=BUILDING, 2=GROUND, 3=TREE).  If configured, wind_speed and
     * wind_velocity are appended as big-endian float cell arrays.
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
    std::vector<float> windSpeedByLevel_;
    std::vector<std::array<float,3> > windVelocityByLevel_;
    std::vector<std::uint32_t> treeWindCellIndices_;
    std::vector<float> treeWindAttenuations_;
    double referenceWindSpeed_ = std::numeric_limits<double>::quiet_NaN();
    double referenceWindDirection_ = std::numeric_limits<double>::quiet_NaN();
    double roughnessLength_ = std::numeric_limits<double>::quiet_NaN();
    double referenceWindHeight_ = 10.;
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
