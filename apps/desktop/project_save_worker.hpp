#pragma once

#include "project_file.hpp"
#include <QMetaObject>
#include <QObject>
#include <QString>
#include <condition_variable>
#include <cstdio>
#include <exception>
#include <functional>
#include <mutex>
#include <optional>
#include <thread>
#include <utility>

namespace hornsculpt {
// One managed save at a time. A persistent worker owns each value snapshot;
// completion is delivered on the owning thread without joining it there.
class ProjectSaveWorker final : public QObject {
  public:
    using Completion = std::function<void(QString)>;
    using Writer = std::function<void(const QString &, const Project &)>;

    explicit ProjectSaveWorker(QObject *parent = nullptr, Writer writer = {})
        : QObject(parent), writer_(writer ? std::move(writer) : Writer{writeProject}) {}

    ~ProjectSaveWorker() override {
        if (thread_.joinable()) {
            thread_.request_stop();
            ready_.notify_all();
            thread_.join();
        }
    }

    bool busy() const noexcept { return busy_; }

    bool start(QString path, Project project, Completion completion) {
        if (busy_)
            return false;
        Task task{std::move(path), std::move(project), std::move(completion)};
        // Construction can throw std::system_error; the caller restores its UI.
        if (!thread_.joinable())
            thread_ = std::jthread([this](std::stop_token stop) { run(stop); });
        {
            std::lock_guard lock(mutex_);
            pending_.emplace(std::move(task));
        }
        busy_ = true;
        ready_.notify_one();
        return true;
    }

  private:
    struct Task {
        QString path;
        Project project;
        Completion completion;
    };

    void run(std::stop_token stop) {
        while (true) {
            std::optional<Task> task;
            {
                std::unique_lock lock(mutex_);
                ready_.wait(lock, stop, [this] { return pending_.has_value(); });
                if (stop.stop_requested() && !pending_)
                    return;
                task = std::move(pending_);
                pending_.reset();
            }
            QString error;
            try {
                writer_(task->path, task->project);
            } catch (const std::exception &e) {
                error = QString::fromUtf8(e.what());
                if (error.isEmpty())
                    error = tr("Fehler beim Speichern.");
            } catch (...) {
                error = tr("Unbekannter Fehler beim Speichern.");
            }
            const bool queued = QMetaObject::invokeMethod(
                this,
                [this, completion = std::move(task->completion), error = std::move(error)]() mutable {
                    busy_ = false;
                    completion(std::move(error));
                },
                Qt::QueuedConnection);
            // A live owner normally always accepts the typed queued functor.
            // Do not silently hide an unexpected delivery failure.
            if (!queued) {
                std::fprintf(stderr, "Project save completion could not be queued.\n");
                return;
            }
        }
    }

    Writer writer_;
    std::mutex mutex_;
    std::condition_variable_any ready_;
    std::optional<Task> pending_;
    std::jthread thread_;
    bool busy_{}; // Accessed only by the owning GUI thread.
};
} // namespace hornsculpt
