#ifndef MASS_GRID_TRANSPORT_H
#define MASS_GRID_TRANSPORT_H

#include <algorithm>
#include <cmath>
#include <vector>

#include "pointer_merged_transport.h"

/**
 * Fixed mass-coordinate water column of one pipe, for the "massGrid" dynamic thermal model.
 *
 * The pipe water is `n` cells of equal mass. They still move with the water, as the parcels of
 * the Lagrangian model do, but their number, mass and storage never change: the column is a ring
 * of n+1 slots and one phase phi in [0,1], the volume (in cells) of the partial slot at the tail.
 * Seen from the tail, the slots hold phi, then n-1 full cells, then 1-phi at the head, so the
 * column always tiles the pipe exactly.
 *
 * Moving the water by a signed displacement r (in cells, r = m_dot dt / dm) is an index rotation
 * plus an update of the boundary slots: floor(phi + r) slots leave at the downstream end and the
 * same number enter at the upstream end at the inlet temperature. Cells are never interpolated
 * into each other, so the transport adds no numerical diffusion; the only averaging is the mixing
 * of the inflow into the partial inlet cell. A reversal only mirrors the direction the index
 * moves: the storage, the cell masses and the wall segments are untouched.
 *
 * The wall segments of the pipe are the n fixed cells of the grid. With the column offset by phi,
 * segment k holds exactly slot k (volume phi) and slot k+1 (volume 1-phi), for every k. That is
 * the fixed two-diagonal "sparse structure" between the moving water and the fixed wall; only the
 * weights phi and 1-phi change with the flow.
 */
class MassGridColumn {
public:
    MassGridColumn() : n(0u), start(0u), phi(0.), mergeBuildCount(0u),
                       mergeGroupCountSum(0u), mergeMaxGroupCount(0u) {}

    /** Allocates the n+1 slots once, at a uniform temperature. Never called again afterwards. */
    void init(unsigned int cellCount, float temperature)
    {
        n = std::max(1u, cellCount);
        slots.assign(n + 1u, temperature);
        start = 0u;
        phi = 0.;
        mergedView.reserve(n + 1u);
        mergeBuildCount = 0u;
        mergeGroupCountSum = 0u;
        mergeMaxGroupCount = 0u;
    }

    bool ready() const { return n > 0u; }
    unsigned int cells() const { return n; }
    double phase() const { return phi; }

    /** Slot i counted from the tail (connectedNode[0]) end, i in [0, n]. */
    float& fromTail(unsigned int i) { return slots[wrap(start + i)]; }
    float fromTail(unsigned int i) const { return slots[wrap(start + i)]; }

    /** Walks the ring from the tail towards the head without a modulo per step. */
    class Cursor {
    public:
        Cursor(float* at, float* first, float* last) : p(at), begin(first), end(last) {}
        float& operator*() const { return *p; }
        void next() { if (++p == end) p = begin; }
    private:
        float* p;
        float* begin;
        float* end;
    };
    Cursor fromTailCursor(unsigned int i)
    {
        return Cursor(slots.data() + wrap(start + i), slots.data(), slots.data() + slots.size());
    }

    /** Volume of slot i, in cells. */
    double slotVolume(unsigned int i) const { return i == 0u ? phi : (i == n ? 1. - phi : 1.); }

    /** Sum of volume x temperature over the column, in cell units [cell.K]. */
    double enthalpy() const
    {
        // Full cells in storage order, then the two partial slots corrected to their volume.
        double h = 0.;
        for (unsigned int j = 0u; j <= n; ++j) h += slots[j];
        return h - (1. - phi) * fromTail(0u) - phi * fromTail(n);
    }

    double meanTemperature() const { return enthalpy() / static_cast<double>(n); }

    /**
     * Temperature at the downstream end of a still column: the last slot that holds water, on the
     * side the flow last left from.
     */
    float endTemperature(bool towardsHead) const
    {
        const double eps = 1.e-9;
        if (towardsHead) return (1. - phi > eps || n == 0u) ? fromTail(n) : fromTail(n - 1u);
        return (phi > eps || n == 0u) ? fromTail(0u) : fromTail(1u);
    }

    /**
     * Moves the water by `r` cells (signed: > 0 from tail to head) with water at `inletTemperature`
     * entering upstream. Returns the mean temperature of the water that left downstream during the
     * move, weighted by volume, so that sum(out) = enthalpy lost by the column + inflow exactly.
     */
    float advect(double r, float inletTemperature)
    {
        if (r > 0.) return advectDirected<true>(r, inletTemperature);
        if (r < 0.) return advectDirected<false>(-r, inletTemperature);
        return endTemperature(true);
    }

    /**
     * Rebuild a transient cascaded merge view, then advect the unchanged physical plug grid.
     * Complete groups are reduced through their enthalpy-weighted temperature; a group intersected
     * by the outlet cut falls back to its original plugs, preserving the exact conservative result.
     */
    float advectPointerMerged(double r, float inletTemperature, float temperatureThreshold)
    {
        if (r > 0.) return advectDirected<true>(r, inletTemperature, temperatureThreshold);
        if (r < 0.) return advectDirected<false>(-r, inletTemperature, temperatureThreshold);
        buildMergedView<true>(temperatureThreshold);
        recordMergedView();
        return endTemperature(true);
    }

    unsigned int currentMergeGroupCount() const { return mergedView.groups(); }
    uint64_t mergedViewBuildCount() const { return mergeBuildCount; }
    uint64_t mergedGroupCountSum() const { return mergeGroupCountSum; }
    unsigned int maxMergeGroupCount() const { return mergeMaxGroupCount; }
    const PointerMergedView& currentMergedView() const { return mergedView; }

private:
    unsigned int wrap(unsigned int j) const { return j >= n + 1u ? j - (n + 1u) : j; }

    // Slot `i` counted from the inlet of the current flow direction.
    template <bool Forward>
    float& fromInlet(unsigned int i) { return fromTail(Forward ? i : n - i); }

    template <bool Forward>
    unsigned int physicalIndexFromInlet(unsigned int i) const
    {
        const unsigned int logical = Forward ? i : n - i;
        return wrap(start + logical);
    }

    template <bool Forward>
    void buildMergedView(float temperatureThreshold)
    {
        mergedView.begin(temperatureThreshold);
        const double psi = Forward ? phi : 1. - phi;
        // Outlet-to-inlet order.  The first and last entries can be partial; zero-volume entries
        // are omitted, so no artificial group bridges an empty boundary slot.
        for (unsigned int j = 0u; j <= n; ++j) {
            const unsigned int inletIndex = n - j;
            const double volume = j == 0u ? 1. - psi : (j == n ? psi : 1.);
            mergedView.append(static_cast<uint32_t>(physicalIndexFromInlet<Forward>(inletIndex)),
                              fromInlet<Forward>(inletIndex), volume);
        }
        mergedView.end();
    }

    void recordMergedView()
    {
        ++mergeBuildCount;
        mergeGroupCountSum += mergedView.groups();
        mergeMaxGroupCount = std::max(mergeMaxGroupCount, mergedView.groups());
    }

    template <bool Forward>
    float advectDirected(double r, float tin, float mergeThreshold = -1.f)
    {
        // Phase of the inlet slot, as seen from the inlet.
        const double psi = Forward ? phi : 1. - phi;
        const double shifted = psi + r;
        const double kReal = std::floor(shifted);
        const double nd = static_cast<double>(n);

        if (mergeThreshold >= 0.f) {
            buildMergedView<Forward>(mergeThreshold);
            recordMergedView();
        }

        if (kReal >= nd) {
            // The advancement reaches the inlet slot itself (or flushes the whole pipe). Whatever
            // is left of the old inlet slot ends up in the new outlet slot; the rest is inflow.
            const double k = kReal;
            const double psiNew = shifted - k;
            const double remnant = std::max(0., nd - r); // old inlet-slot water still in the pipe
            const float oldInlet = fromInlet<Forward>(0u);
            const double hOut = mergeThreshold >= 0.f
                              ? mergedView.prefixEnthalpy(std::min(r, nd), slots) + std::max(0., r - nd) * tin
                              : enthalpy() - remnant * oldInlet + std::max(0., r - nd) * tin;
            for (unsigned int i = 0u; i <= n; ++i) slots[i] = tin;
            const double outletVolume = 1. - psiNew;
            if (outletVolume > 0.) {
                fromInlet<Forward>(n) = static_cast<float>((remnant * oldInlet + (outletVolume - remnant) * tin) / outletVolume);
            }
            setInletPhase<Forward>(psiNew);
            return static_cast<float>(hOut / r);
        }

        const unsigned int k = static_cast<unsigned int>(kReal);
        const double psiNew = shifted - kReal;

        // Outflow: the outlet slot (1 - psi), then k-1 full cells, then psiNew of the next one.
        double hOut;
        if (mergeThreshold >= 0.f) {
            hOut = mergedView.prefixEnthalpy(r, slots);
        }
        else if (k == 0u) {
            hOut = r * fromInlet<Forward>(n);
        }
        else {
            hOut = (1. - psi) * fromInlet<Forward>(n);
            for (unsigned int j = 1u; j < k; ++j) hOut += fromInlet<Forward>(n - j);
            hOut += psiNew * fromInlet<Forward>(n - k);
        }

        // Inflow: fill the inlet slot, then k-1 full new cells and a new partial inlet slot.
        float& inlet = fromInlet<Forward>(0u);
        if (k == 0u) {
            inlet = static_cast<float>((psi * inlet + r * tin) / shifted);
        }
        else {
            inlet = static_cast<float>(psi * inlet + (1. - psi) * tin);
            // Rotate: the k drained slots become the k new slots at the inlet.
            start = Forward ? wrap(start + (n + 1u) - k) : wrap(start + k);
            for (unsigned int j = 0u; j < k; ++j) fromInlet<Forward>(j) = tin;
        }
        setInletPhase<Forward>(psiNew);
        return static_cast<float>(hOut / r);
    }

    template <bool Forward>
    void setInletPhase(double psiNew) { phi = Forward ? psiNew : 1. - psiNew; }

    std::vector<float> slots;
    unsigned int n;
    unsigned int start;
    double phi;
    PointerMergedView mergedView;
    uint64_t mergeBuildCount;
    uint64_t mergeGroupCountSum;
    unsigned int mergeMaxGroupCount;
};

#endif
