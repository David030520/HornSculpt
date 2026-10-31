#pragma once

#include <cstddef>
#include <hornsculpt/draft_geometry.hpp>
#include <stop_token>
#include <string>

namespace hornsculpt {

enum class DraftStepExportKind { solid, innerSurface, sections };

struct DraftStepExportOptions {
    DraftStepExportKind kind{DraftStepExportKind::solid};
    double toleranceMetres{0.00005};
    // Includes entry and mouth. Sections also exports all four cardinal profiles.
    std::size_t sectionCount{9};
};

struct DraftStepExportResult {
    std::string data;
    // Maximum measured 3D discrepancy on independent parameter-cell samples;
    // this is not a certified global Hausdorff bound.
    double maximumDeviationMetres{};
    std::size_t surfaceCount{}, curveCount{};
};

// STEP uses explicit metre units. The material solid spans z=0..L, with the
// Draft's radial wall thickness and annular ends: both acoustic openings stay
// open. No driven disk or simulation rear closure is exported. Surfaces and
// curves are interpolated B-splines, never converted preview triangles.
// Throws on invalid topology, failed accuracy checks, cancellation or budget.
DraftStepExportResult exportDraftStep(const DraftGeometry &geometry,
                                      const DraftStepExportOptions &options = {},
                                      std::stop_token stop = {});

} // namespace hornsculpt
