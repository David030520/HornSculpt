#include "geometry_mesh_worker.hpp"
#include <QCoreApplication>
#include <QElapsedTimer>
#include <QEventLoop>
#include <QProcess>
#include <QSemaphore>
#include <QThread>
#include <QTimer>
#ifdef Q_OS_WIN
#define NOMINMAX
#include <windows.h>
#endif
#include <algorithm>
#include <atomic>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <stop_token>

using hornsculpt::GeometryMeshWorker;
void require(bool ok, const char *message) {
    if (!ok)
        throw std::runtime_error(message);
}
void waitFor(const std::function<bool()> &condition) {
    QElapsedTimer t;
    t.start();
    while (!condition() && t.elapsed() < 10000) {
        QCoreApplication::processEvents(QEventLoop::AllEvents, 5);
        QThread::yieldCurrentThread();
    }
    require(condition(), "Worker contract timeout");
}
int main(int argc, char **argv) {
    QCoreApplication app(argc, argv);
    try {
        if (app.arguments().contains("--reentrant-child")) {
            GeometryMeshWorker child;
            QSemaphore entered, release;
            child.submit([&](std::stop_token token) {
                std::stop_callback reentrant(token, [&] { child.request_stop(); });
                entered.release();
                release.acquire();
            });
            waitFor([&] { return entered.available() > 0; });
            child.request_stop();
            release.release();
            child.shutdown();
            return 0;
        }
        QProcess child;
#ifdef Q_OS_WIN
        child.setCreateProcessArgumentsModifier(
            [](QProcess::CreateProcessArguments *args) { args->flags |= CREATE_NO_WINDOW; });
#endif
        child.start(QCoreApplication::applicationFilePath(), {"--reentrant-child"});
        require(child.waitForStarted(5000), "Cannot start reentrancy test");
        if (!child.waitForFinished(5000)) {
            child.kill();
            child.waitForFinished(5000);
            throw std::runtime_error("Reentrant stop callback deadlocked");
        }
        require(child.exitStatus() == QProcess::NormalExit && child.exitCode() == 0,
                "Reentrant stop callback failed");
        GeometryMeshWorker worker;
        QSemaphore entered, release, latestFinished;
        std::atomic<int> active{}, maxActive{}, ranMiddle{}, ranLatest{}, oldPublished{};
        worker.submit([&](std::stop_token stop) {
            const int n = ++active;
            maxActive.store(std::max(maxActive.load(), n));
            entered.release();
            release.acquire();
            if (!stop.stop_requested())
                ++oldPublished;
            --active;
        });
        waitFor([&] { return entered.available() > 0; });
        worker.submit([&](std::stop_token) { ++ranMiddle; });
        worker.submit([&](std::stop_token) {
            const int n = ++active;
            maxActive.store(std::max(maxActive.load(), n));
            ++ranLatest;
            --active;
            latestFinished.release();
        });
        bool eventLoopAlive = false;
        QTimer::singleShot(0, [&] { eventLoopAlive = true; });
        QCoreApplication::processEvents();
        const bool pendingDidNotStart = ranMiddle.load() == 0 && ranLatest.load() == 0;
        release.release();
        waitFor([&] { return latestFinished.available() > 0; });
        require(eventLoopAlive && pendingDidNotStart, "GUI did not return while old job was gated");
        require(oldPublished == 0 && ranMiddle == 0 && ranLatest == 1 && maxActive == 1,
                "Cancelled job published, pending not replaced, or jobs overlapped");
        QSemaphore entered2, release2, done2, pendingDone;
        std::atomic<int> pendingRan{}, pendingStopped{}, stopped{};
        worker.submit([&](std::stop_token token) {
            entered2.release();
            release2.acquire();
            stopped = token.stop_requested();
            done2.release();
        });
        waitFor([&] { return entered2.available() > 0; });
        worker.submit([&](std::stop_token token) {
            ++pendingRan;
            pendingStopped = token.stop_requested();
            pendingDone.release();
        });
        worker.request_stop();
        release2.release();
        waitFor([&] { return done2.available() > 0 && pendingDone.available() > 0; });
        worker.shutdown();
        require(stopped == 1 && pendingRan == 1 && pendingStopped == 1,
                "Cancel did not stop active and pending completion paths");
        bool rejected = false;
        try {
            worker.submit([](std::stop_token) {});
        } catch (const std::logic_error &) {
            rejected = true;
        }
        require(rejected, "Submit after shutdown was not rejected");
        GeometryMeshWorker shuttingDown;
        QSemaphore shutdownEntered, shutdownStopped, shutdownRelease;
        std::atomic<bool> shutdownFinished{}, shutdownPublished{};
        shuttingDown.submit([&](std::stop_token token) {
            std::stop_callback onStop(token, [&] { shutdownStopped.release(); });
            shutdownEntered.release();
            shutdownRelease.acquire();
            if (!token.stop_requested())
                shutdownPublished = true;
            shutdownFinished = true;
        });
        waitFor([&] { return shutdownEntered.available() > 0; });
        std::jthread releaser([&] {
            shutdownStopped.acquire();
            shutdownRelease.release();
        });
        shuttingDown.shutdown();
        releaser.join();
        require(shutdownFinished && !shutdownPublished,
                "Shutdown returned early or allowed stale publication");
        std::cout << "Geometry/mesh worker cancellation and latest-pending contract passed.\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
