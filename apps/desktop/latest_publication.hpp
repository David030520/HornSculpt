#pragma once

#include <mutex>
#include <optional>
#include <utility>

namespace hornsculpt {
// One queued delivery at most; later values replace an undelivered value.
// The owner must arrange a delivery whenever offer() returns true.
template <typename T> class LatestPublication {
  public:
    bool offer(T value) {
        std::lock_guard lock(mutex_);
        value_ = std::move(value);
        if (pending_) return false;
        pending_ = true;
        return true;
    }
    std::optional<T> take() {
        std::lock_guard lock(mutex_);
        pending_ = false;
        return std::exchange(value_, {});
    }
  private:
    std::mutex mutex_;
    std::optional<T> value_;
    bool pending_{};
};
} // namespace hornsculpt
