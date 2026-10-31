#pragma once

#include <hornsculpt/lem.hpp>
#include <numbers>

namespace hornsculpt {

// Shared editable concept seed, not a verified dry-mass driver model.
// VISATON TI 100 publishes Mms=4 g, fs=62 Hz, Qms=5.95; Cms and Rms are derived
// from that consistent triplet. Mmd=4 g is explicitly an assumption, without
// an unverified radiation-mass correction. UI/model names must identify it.
// https://www.visaton.de/de/produkte/chassis/tiefmitteltoener/ti-100-8-ohm
inline lem::DriverParameters ti100ConceptParameters() {
    const double resonance = 2 * std::numbers::pi * 62;
    return {6.4,   .0007, 4.98, .004, 1 / (resonance * resonance * .004), resonance * .004 / 5.95,
            .00544};
}

} // namespace hornsculpt
