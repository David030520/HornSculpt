#pragma once

#include "project_save_worker.hpp"
#include <QDateTime>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonDocument>
#include <QJsonObject>
#include <QLockFile>
#include <QSaveFile>
#include <QUuid>
#include <memory>
#include <stdexcept>
#include <vector>

namespace hornsculpt {
// Each process owns unique recovery files and an exclusive session lock.
// Normal project files are never written through this store.
class RecoveryStore final {
  public:
    struct Entry {
        QString path;
        QString originalPath;
        QDateTime modified;
        std::unique_ptr<QLockFile> lock;
    };

    explicit RecoveryStore(QString directory, ProjectSaveWorker::Writer writer = {})
        : directory_(std::move(directory)),
          worker_(std::make_unique<ProjectSaveWorker>(nullptr, std::move(writer))) {
        if (!QDir().mkpath(directory_))
            throw std::runtime_error("Wiederherstellungsordner konnte nicht angelegt werden.");
        newSession();
    }
    ~RecoveryStore() {
        worker_.reset(); // On orderly shutdown, finish an in-flight atomic write before cleanup.
        for (const auto &record : records_)
            if (record->retired) remove(*record);
    }
    RecoveryStore(const RecoveryStore &) = delete;
    RecoveryStore &operator=(const RecoveryStore &) = delete;

    bool busy() const { return worker_->busy(); }
    void newSession() {
        retire();
        auto record = std::make_shared<Record>();
        record->path = QDir(directory_).filePath(QUuid::createUuid().toString(QUuid::WithoutBraces) + ".hsproj");
        record->lock = std::make_unique<QLockFile>(record->path + ".lock");
        record->lock->setStaleLockTime(0);
        if (!record->lock->tryLock(0))
            throw std::runtime_error("Sitzungssperre für die Wiederherstellung fehlgeschlagen.");
        records_.push_back(record);
        current_ = std::move(record);
    }
    void retire() {
        if (!current_) return;
        current_->retired = true;
        if (!busy()) remove(*current_);
        current_.reset();
    }
    bool save(Project project, const QString &originalPath, std::function<void(QString)> done) {
        if (busy() || !current_) return false;
        auto record = current_;
        if (!record->metadataWritten) {
            QSaveFile meta(record->path + ".json");
            if (!meta.open(QIODevice::WriteOnly))
                throw std::runtime_error("Wiederherstellungsmetadaten konnten nicht geschrieben werden.");
            const auto bytes = QJsonDocument(QJsonObject{{"originalPath", originalPath}}).toJson();
            if (meta.write(bytes) != bytes.size() || !meta.commit())
                throw std::runtime_error("Wiederherstellungsmetadaten konnten nicht gespeichert werden.");
            record->metadataWritten = true;
        }
        return worker_->start(record->path, std::move(project),
            [record, done = std::move(done)](QString error) mutable {
                if (record->retired) remove(*record);
                if (done) done(std::move(error));
            });
    }
    static std::vector<Entry> available(const QString &directory) {
        std::vector<Entry> result;
        const QDir dir(directory);
        for (const auto &info : dir.entryInfoList({"*.hsproj"}, QDir::Files, QDir::Time)) {
            auto lock = std::make_unique<QLockFile>(info.absoluteFilePath() + ".lock");
            lock->setStaleLockTime(0);
            if (!lock->tryLock(0)) continue;
            QString original;
            QFile meta(info.absoluteFilePath() + ".json");
            if (meta.open(QIODevice::ReadOnly))
                original = QJsonDocument::fromJson(meta.readAll()).object().value("originalPath").toString();
            result.push_back({info.absoluteFilePath(), original, info.lastModified(), std::move(lock)});
        }
        return result;
    }
    static void discard(Entry &entry) {
        QFile::remove(entry.path);
        QFile::remove(entry.path + ".json");
        entry.lock.reset();
    }

  private:
    struct Record {
        QString path;
        std::unique_ptr<QLockFile> lock;
        bool retired{}, metadataWritten{};
    };
    static void remove(Record &record) {
        QFile::remove(record.path);
        QFile::remove(record.path + ".json");
        record.lock.reset();
    }
    QString directory_;
    std::unique_ptr<ProjectSaveWorker> worker_;
    std::shared_ptr<Record> current_;
    std::vector<std::shared_ptr<Record>> records_;
};
} // namespace hornsculpt
