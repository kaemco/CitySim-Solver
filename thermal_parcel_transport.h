#ifndef THERMAL_PARCEL_TRANSPORT_H
#define THERMAL_PARCEL_TRANSPORT_H

#include <algorithm>
#include <cmath>
#include <vector>

// Parcel handling shared by the dynamic thermal transport models.

/** Parcels closer than this in temperature may share one parcel [K]. */
static const float THERMAL_PARCEL_MERGE_TOLERANCE = 0.01f;

// Merges `volume` of water at `temperature` into `into` when within `tolerance` and `volumeCap` (enthalpy-weighted).
template <class Parcel>
inline bool thermalMergeParcel(Parcel& into, float volume, float temperature,
                               float tolerance, float volumeCap)
{
    if (std::fabs(into.temperature - temperature) >= tolerance) return false;
    const double mergedVolume = static_cast<double>(into.volume) + volume;
    if (mergedVolume > volumeCap * (1. + 1.e-6)) return false;
    if (mergedVolume > 0.) {
        into.temperature = static_cast<float>((static_cast<double>(into.volume) * into.temperature
                                               + static_cast<double>(volume) * temperature) / mergedVolume);
    }
    into.volume = static_cast<float>(mergedVolume);
    return true;
}

/**
 * Walks the ordered parcels of a pipe along its fixed wall segments, starting from the tail end.
 *
 * The segment is tracked by index and remaining room, not recomputed as position/segmentVolume:
 * on an exact boundary that division can truncate to the previous segment, which has no room
 * left, and the whole remainder of a parcel then landed in that one segment. Any rounding excess
 * of the water column goes to the last segment.
 */
class ThermalSegmentCursor {
public:
    ThermalSegmentCursor(unsigned int segmentCount, float segmentVolume) :
        count(std::max(1u, segmentCount)), volume(segmentVolume), index(0u), room(segmentVolume) {}

    /** Volume of the next slice of a parcel that still has `remaining` to place, lying in segment(). */
    float nextSlice(float remaining)
    {
        if (index + 1u < count && room <= volume * 1.e-6f) {
            ++index;
            room = volume;
        }
        if (index + 1u >= count) return remaining;
        return std::min(remaining, room);
    }

    void consume(float slice) { room -= slice; }
    unsigned int segment() const { return index; }

private:
    unsigned int count;
    float volume;
    unsigned int index;
    float room;
};

// Rebuilds the parcel sequence of a pipe from the heated slices: split at `tolerance`, merge up to `mergeVolumeCap`.
template <class Parcel>
class ThermalSliceRegrouper {
public:
    ThermalSliceRegrouper(std::vector<Parcel>& output, float mergeTolerance, float volumeCap, bool keepParcels) :
        out(output), tolerance(mergeTolerance), mergeVolumeCap(volumeCap), preserveParcels(keepParcels),
        source(), parcelVolume(0.), pieceVolume(0.), pieceEnthalpy(0.), emittedVolume(0.), firstPiece(true) {}

    void beginParcel(const Parcel& parcel)
    {
        source = parcel;
        parcelVolume = parcel.volume;
        pieceVolume = 0.;
        pieceEnthalpy = 0.;
        emittedVolume = 0.;
        firstPiece = true;
    }

    void addSlice(float volume, float temperature)
    {
        if (!preserveParcels && pieceVolume > 0. &&
            std::fabs(temperature - pieceEnthalpy / pieceVolume) >= tolerance) {
            emit(pieceVolume, pieceEnthalpy / pieceVolume);
            pieceVolume = 0.;
            pieceEnthalpy = 0.;
        }
        pieceVolume += volume;
        pieceEnthalpy += static_cast<double>(volume) * temperature;
    }

    void endParcel()
    {
        if (pieceVolume <= 0.) return;
        // The last piece takes the exact rest of the parcel, so that splitting a parcel never
        // changes the volume of the column by rounding.
        double volume = parcelVolume - emittedVolume;
        if (volume <= 0.) volume = pieceVolume;
        emit(volume, pieceEnthalpy / pieceVolume);
    }

private:
    void emit(double volume, double temperature)
    {
        const bool mayJoinPrevious = !preserveParcels && firstPiece && !out.empty();
        if (!(mayJoinPrevious && thermalMergeParcel(out.back(), static_cast<float>(volume),
                                                    static_cast<float>(temperature), tolerance, mergeVolumeCap))) {
            Parcel piece = source;
            piece.temperature = static_cast<float>(temperature);
            piece.volume = static_cast<float>(volume);
            out.push_back(piece);
        }
        emittedVolume += volume;
        firstPiece = false;
    }

    std::vector<Parcel>& out;
    float tolerance;
    float mergeVolumeCap;
    bool preserveParcels;
    Parcel source;
    double parcelVolume;
    double pieceVolume;
    double pieceEnthalpy;
    double emittedVolume;
    bool firstPiece;
};

#endif
