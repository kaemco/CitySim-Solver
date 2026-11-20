#ifndef POINTER_MERGED_TRANSPORT_H
#define POINTER_MERGED_TRANSPORT_H

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

/**
 * A transient CRS view of temperature-compatible, adjacent Lagrangian plugs.
 *
 * The physical plugs live elsewhere and are never merged.  `plugIndices` is the ordered list of
 * stable integer references into that storage and `groupOffsets` delimits consecutive groups in
 * the list.  The remaining arrays are structure-of-arrays group data.  This is the pointer-list /
 * merged-pointer-list requested by the algorithm, expressed as indices rather than C++ pointers
 * so storage can be copied to an accelerator without pointer relocation.
 *
 * Groups are rebuilt from scratch before every advection.  A candidate joins the current group
 * when its temperature differs from the current enthalpy-weighted group mean by at most epsilon.
 * Updating the mean after every accepted plug makes the merge cascaded and deterministic.
 */
class PointerMergedView {
public:
    PointerMergedView() : open(false), epsilon(0.f) {}

    void reserve(unsigned int plugCapacity)
    {
        plugIndices.reserve(plugCapacity);
        plugVolumes.reserve(plugCapacity);
        groupOffsets.reserve(plugCapacity + 1u);
        groupVolumes.reserve(plugCapacity);
        groupTemperatures.reserve(plugCapacity);
    }

    void begin(float temperatureThreshold)
    {
        plugIndices.clear();
        plugVolumes.clear();
        groupOffsets.clear();
        groupVolumes.clear();
        groupTemperatures.clear();
        open = false;
        epsilon = std::max(0.f, temperatureThreshold);
    }

    /** Add one positive-volume plug in advection order (outlet towards inlet). */
    void append(uint32_t physicalIndex, float temperature, double volume)
    {
        if (volume <= 0.) return;

        if (!open || std::fabs(static_cast<double>(temperature) - groupTemperatures.back()) > epsilon) {
            groupOffsets.push_back(static_cast<uint32_t>(plugIndices.size()));
            groupVolumes.push_back(0.);
            groupTemperatures.push_back(0.);
            open = true;
        }

        plugIndices.push_back(physicalIndex);
        plugVolumes.push_back(volume);

        const size_t group = groupVolumes.size() - 1u;
        const double oldVolume = groupVolumes[group];
        const double newVolume = oldVolume + volume;
        groupTemperatures[group] = (groupTemperatures[group] * oldVolume
                                  + static_cast<double>(temperature) * volume) / newVolume;
        groupVolumes[group] = newVolume;
    }

    void end()
    {
        if (open) groupOffsets.push_back(static_cast<uint32_t>(plugIndices.size()));
        open = false;
    }

    unsigned int plugs() const { return static_cast<unsigned int>(plugIndices.size()); }
    unsigned int groups() const { return static_cast<unsigned int>(groupVolumes.size()); }

    const std::vector<uint32_t>& indices() const { return plugIndices; }
    const std::vector<uint32_t>& offsets() const { return groupOffsets; }
    const std::vector<double>& volumes() const { return groupVolumes; }
    const std::vector<double>& temperatures() const { return groupTemperatures; }

    /**
     * Integrate volume*temperature over a prefix in advection order.  Complete CRS rows use their
     * enthalpy-weighted mean in O(1); only a row cut by the requested volume visits its plugs.
     * The physical plug values are untouched, including in the cut-row path.
     */
    double prefixEnthalpy(double requestedVolume,
                          const std::vector<float>& physicalTemperatures) const
    {
        double remaining = std::max(0., requestedVolume);
        double enthalpy = 0.;
        const double tolerance = 1.e-12 * std::max(1., requestedVolume);

        for (size_t group = 0u; group < groupVolumes.size() && remaining > tolerance; ++group) {
            const double groupVolume = groupVolumes[group];
            if (remaining + tolerance >= groupVolume) {
                enthalpy += groupVolume * groupTemperatures[group];
                remaining -= groupVolume;
                continue;
            }

            const uint32_t first = groupOffsets[group];
            const uint32_t last = groupOffsets[group + 1u];
            for (uint32_t member = first; member < last && remaining > tolerance; ++member) {
                const double take = std::min(remaining, plugVolumes[member]);
                enthalpy += take * physicalTemperatures[plugIndices[member]];
                remaining -= take;
            }
            break;
        }
        return enthalpy;
    }

private:
    // CRS connectivity: group g owns plugIndices[groupOffsets[g]..groupOffsets[g+1]).
    std::vector<uint32_t> plugIndices;
    std::vector<uint32_t> groupOffsets;

    // Slot volumes accompany the integer pointers so a cut row remains exact at partial slots.
    std::vector<double> plugVolumes;

    // Structure-of-arrays aggregate state, suitable for a segmented scan/reduction on a GPU.
    std::vector<double> groupVolumes;
    std::vector<double> groupTemperatures;

    bool open;
    float epsilon;
};

#endif
