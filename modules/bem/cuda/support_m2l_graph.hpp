#pragma once

#include <chrono>
#include <cstddef>
#include <cuda_runtime_api.h>
#include <stdexcept>
#include <string>
#include <utility>

namespace hornsculpt::wideband::detail {
// Borrowed-stream graph owner. Callers drain the stream before invalidating
// any captured pointer. Capture does not execute or count numerical work.
class M2lChunkGraph {
  public:
    struct Preparation {
        std::size_t nodes{};
        double captureMilliseconds{}, instantiationMilliseconds{};
    };
    M2lChunkGraph() = default;
    ~M2lChunkGraph() {
        if (executable_)
            cudaGraphExecDestroy(executable_);
    }
    M2lChunkGraph(const M2lChunkGraph &) = delete;
    M2lChunkGraph &operator=(const M2lChunkGraph &) = delete;
    M2lChunkGraph(M2lChunkGraph &&other) noexcept
        : executable_(std::exchange(other.executable_, nullptr)) {}
    M2lChunkGraph &operator=(M2lChunkGraph &&other) noexcept {
        if (this != &other) {
            if (executable_)
                cudaGraphExecDestroy(executable_);
            executable_ = std::exchange(other.executable_, nullptr);
        }
        return *this;
    }
    explicit operator bool() const noexcept {
        return executable_ != nullptr;
    }
    template <class Enqueue>
    bool prepare(cudaStream_t stream, std::size_t maximumNodes, Enqueue enqueue,
                 Preparation &timing) {
        if (executable_)
            throw std::logic_error("M2L graph is already instantiated.");
        struct Captured {
            cudaGraph_t value{};
            ~Captured() {
                if (value)
                    cudaGraphDestroy(value);
            }
        } graph;
        const auto begin = std::chrono::steady_clock::now();
        check(cudaStreamBeginCapture(stream, cudaStreamCaptureModeThreadLocal),
              "M2L graph capture begin");
        try {
            enqueue();
        } catch (...) {
            // End even an invalidated capture before the outer operation drains
            // its stream. Cancellation cannot leave the borrowed stream capturing.
            const auto status = cudaStreamEndCapture(stream, &graph.value);
            if (status != cudaSuccess && status != cudaErrorStreamCaptureInvalidated)
                check(status, "M2L graph capture failure cleanup");
            throw;
        }
        check(cudaStreamEndCapture(stream, &graph.value), "M2L graph capture end");
        check(cudaGraphGetNodes(graph.value, nullptr, &timing.nodes), "M2L graph node count");
        timing.captureMilliseconds = elapsed(begin);
        if (timing.nodes > maximumNodes)
            return false;
        const auto instantiate = std::chrono::steady_clock::now();
        check(cudaGraphInstantiateWithFlags(&executable_, graph.value, 0),
              "M2L graph instantiation");
        timing.instantiationMilliseconds = elapsed(instantiate);
        return true;
    }
    void launch(cudaStream_t stream) const {
        if (!executable_)
            throw std::logic_error("M2L graph launch without an executable.");
        check(cudaGraphLaunch(executable_, stream), "M2L graph launch");
    }

  private:
    static void check(cudaError_t status, const char *operation) {
        if (status != cudaSuccess)
            throw std::runtime_error(std::string(operation) + ": " + cudaGetErrorString(status));
    }
    static double elapsed(std::chrono::steady_clock::time_point begin) {
        return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - begin)
            .count();
    }
    cudaGraphExec_t executable_{};
};
} // namespace hornsculpt::wideband::detail