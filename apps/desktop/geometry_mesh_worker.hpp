#pragma once

#include <condition_variable>
#include <functional>
#include <mutex>
#include <stdexcept>
#include <stop_token>
#include <thread>
#include <utility>

namespace hornsculpt {

// Serializes preview or meshing jobs without joining an obsolete job on the UI thread.
// A newer snapshot replaces any waiting snapshot and cancels the active one.
// Jobs and stop callbacks must not let exceptions escape; submit and shutdown
// are called by the owning UI thread, and shutdown joins only at destruction.
class GeometryMeshWorker final {
  public:
    using Job = std::function<void(std::stop_token)>;

    GeometryMeshWorker() = default;
    GeometryMeshWorker(const GeometryMeshWorker &) = delete;
    GeometryMeshWorker &operator=(const GeometryMeshWorker &) = delete;
    ~GeometryMeshWorker() { shutdown(); }

    void submit(Job job) {
        std::stop_source nextStop;
        std::stop_source active{std::nostopstate};
        std::stop_source pending{std::nostopstate};
        {
            std::lock_guard lock(mutex_);
            if (stopping_)
                throw std::logic_error("Cannot submit a job after worker shutdown");
            if (!thread_.joinable())
                thread_ = std::jthread([this] { run(); });
            active = activeStop_;
            pending = pendingStop_;
            pendingStop_ = std::move(nextStop);
            pending_ = std::move(job);
        }
        active.request_stop();
        pending.request_stop();
        ready_.notify_one();
    }

    void request_stop() {
        std::stop_source active{std::nostopstate};
        std::stop_source pending{std::nostopstate};
        {
            std::lock_guard lock(mutex_);
            active = activeStop_;
            pending = pendingStop_;
        }
        active.request_stop();
        pending.request_stop();
    }

    void shutdown() {
        std::stop_source active{std::nostopstate};
        std::stop_source pending{std::nostopstate};
        {
            std::lock_guard lock(mutex_);
            stopping_ = true;
            active = activeStop_;
            pending = pendingStop_;
            pending_ = {};
        }
        active.request_stop();
        pending.request_stop();
        ready_.notify_one();
        if (thread_.joinable())
            thread_.join();
    }

  private:
    void run() {
        for (;;) {
            Job job;
            std::stop_token token;
            {
                std::unique_lock lock(mutex_);
                ready_.wait(lock, [this] { return stopping_ || bool(pending_); });
                if (stopping_)
                    return;
                job = std::move(pending_);
                pending_ = {};
                activeStop_ = pendingStop_;
                token = activeStop_.get_token();
            }
            job(token);
        }
    }

    std::mutex mutex_;
    std::condition_variable ready_;
    Job pending_;
    std::stop_source activeStop_;
    std::stop_source pendingStop_;
    bool stopping_{};
    std::jthread thread_;
};

} // namespace hornsculpt
