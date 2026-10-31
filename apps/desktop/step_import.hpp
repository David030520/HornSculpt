#pragma once
#include <hornsculpt/freecad_link.hpp>
#include <stop_token>

namespace hornsculpt {
std::shared_ptr<const CadSnapshot> readStepSnapshot(const QString &path, double tolerance,
                                                    std::stop_token stop = {});
}
