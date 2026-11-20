#ifndef THERMAL_PARCEL_TRANSPORT_H
#define THERMAL_PARCEL_TRANSPORT_H

#include <algorithm>
#include <vector>

// Parcel handling shared by the dynamic thermal transport models.

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

// Gives a parcel the volume-weighted mean temperature of its heated slices (no merge, no split).
template <class Parcel>
class ThermalParcelAverager {
public:
    explicit ThermalParcelAverager(std::vector<Parcel>& output) :
        out(output), source(), parcelVolume(0.), pieceVolume(0.), pieceEnthalpy(0.) {}

    void beginParcel(const Parcel& parcel)
    {
        source = parcel;
        parcelVolume = parcel.volume;
        pieceVolume = 0.;
        pieceEnthalpy = 0.;
    }

    void addSlice(float volume, float temperature)
    {
        pieceVolume += volume;
        pieceEnthalpy += static_cast<double>(volume) * temperature;
    }

    void endParcel()
    {
        if (pieceVolume <= 0.) return;
        Parcel parcel = source;
        parcel.temperature = static_cast<float>(pieceEnthalpy / pieceVolume);
        parcel.volume = static_cast<float>(parcelVolume);
        out.push_back(parcel);
    }

private:
    std::vector<Parcel>& out;
    Parcel source;
    double parcelVolume;
    double pieceVolume;
    double pieceEnthalpy;
};

#endif
