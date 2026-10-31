#pragma once

namespace hornsculpt {
struct DraftParameters;
// Validate the physical inlet composition before constructing profile providers.
void validateDiffractionSlot(const DraftParameters &parameters);
} // namespace hornsculpt
