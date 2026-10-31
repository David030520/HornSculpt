#include "recovery_store.hpp"
#include <QCoreApplication>
#include <QElapsedTimer>
#include <QCryptographicHash>
#include <QSemaphore>
#include <QTemporaryDir>
#include <QThread>
#include <cstdio>
#include <stdexcept>

using namespace hornsculpt;
namespace {
void require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}
void waitUntil(const std::function<bool()> &condition) {
    QElapsedTimer timer;
    timer.start();
    while (!condition()) {
        if (timer.elapsed() > 10000) throw std::runtime_error("Recovery worker timed out");
        QCoreApplication::processEvents();
        QThread::msleep(2);
    }
}
}
int main(int argc, char **argv) {
    QCoreApplication app(argc, argv);
    try {
        QTemporaryDir directory;
        require(directory.isValid(), "Recovery fixture directory unavailable");
        QString completedError;
        {
            RecoveryStore store(directory.path());
            require(store.save(Project{}, "original.hsproj", [&](QString error) { completedError = error; }),
                    "Recovery write did not start");
            waitUntil([&] { return !store.busy(); });
            require(completedError.isEmpty(), "Recovery write failed");
            require(QDir(directory.path()).entryList({"*.lock"}, QDir::Files).size() == 1,
                    "Recovery session lock missing");
            // Windows locks the lockfile against timestamp edits. Exercise an
            // actually live session beyond QLockFile's default stale interval.
            QThread::sleep(31);
            require(RecoveryStore::available(directory.path()).empty(),
                    "Another instance can claim an old but live recovery session");
        }
        auto recovered = RecoveryStore::available(directory.path());
        require(recovered.size() == 1 && recovered.front().originalPath == "original.hsproj" &&
                    recovered.front().modified.isValid(), "Abandoned recovery was not discoverable");
        (void)readProject(recovered.front().path);
        RecoveryStore::discard(recovered.front());
        require(RecoveryStore::available(directory.path()).empty(), "Discard retained the recovery");

        {
            unsigned writes = 0;
            RecoveryStore store(directory.path(), [&](const QString &path, const Project &project) {
                if (++writes == 2) throw std::runtime_error("Expected disk failure");
                writeProject(path, project);
            });
            QString error;
            require(store.save(Project{}, "original.hsproj", [&](QString value) { error = value; }),
                    "Initial recovery did not start");
            waitUntil([&] { return !store.busy(); });
            require(error.isEmpty(), "Initial recovery failed");
            const auto files = QDir(directory.path()).entryInfoList({"*.hsproj"}, QDir::Files);
            require(files.size() == 1, "Initial recovery is missing");
            QFile file(files.front().absoluteFilePath());
            require(file.open(QIODevice::ReadOnly), "Cannot read previous recovery");
            const auto before = QCryptographicHash::hash(file.readAll(), QCryptographicHash::Sha256);
            file.close();
            require(store.save(Project{}, "original.hsproj", [&](QString value) { error = value; }),
                    "Failing recovery did not start");
            waitUntil([&] { return !store.busy(); });
            require(error.contains("Expected disk failure"), "Recovery failure was hidden");
            require(file.open(QIODevice::ReadOnly) &&
                        QCryptographicHash::hash(file.readAll(), QCryptographicHash::Sha256) == before,
                    "Failed recovery destroyed the previous valid snapshot");
            file.close();
            store.retire();
        }

        QSemaphore entered, release;
        {
            RecoveryStore store(directory.path(), [&](const QString &path, const Project &project) {
                entered.release();
                release.acquire();
                writeProject(path, project);
            });
            struct ReleaseOnExit { QSemaphore &semaphore; ~ReleaseOnExit() { semaphore.release(); } } guard{release};
            require(store.save(Project{}, "first.hsproj", {}), "Delayed recovery did not start");
            require(entered.tryAcquire(1, 5000), "Recovery writer did not enter");
            store.newSession(); // Save/discard/project switch while an old write is in flight.
            release.release();
            waitUntil([&] { return !store.busy(); });
            require(RecoveryStore::available(directory.path()).empty(),
                    "Late completion resurrected a retired recovery");
        }
        std::puts("Recovery sessions, atomic project data, isolation and late-write cleanup passed.");
        return 0;
    } catch (const std::exception &error) {
        std::fprintf(stderr, "%s\n", error.what());
        return 1;
    }
}
