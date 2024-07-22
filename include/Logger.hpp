#pragma once

// Standard
#include <chrono>
#include <iomanip>
#include <iostream>
#include <map>
#include <mutex>
#include <string>

// seqan3
#include <seqan3/core/debug_stream.hpp>

enum class LogLevel { DEBUG, INFO, WARNING, ERROR };

// TODO Change implementation to use os_syncstream
class Logger {
   public:
    static Logger &getInstance() {
        static Logger instance;  // Singleton instance
        return instance;
    }

    static void setLogLevel(const std::string &logLevelString) {
        static const std::map<std::string, LogLevel> stringToLogLevelMap{
            {"debug", LogLevel::DEBUG},
            {"info", LogLevel::INFO},
            {"warning", LogLevel::WARNING},
            {"error", LogLevel::ERROR}};

        auto it = stringToLogLevelMap.find(logLevelString);
        if (it != stringToLogLevelMap.end()) {
            getInstance().logLevel = it->second;
        } else {
            log(LogLevel::ERROR, "Invalid log level: ", logLevelString);
        }
    }

    static void setLogLevel(LogLevel level) { getInstance().logLevel = level; }

    template <typename... Args>
    static void log(LogLevel level, Args &&...args) {
        std::lock_guard<std::mutex> lock(getInstance().logMutex);
        if (level >= getInstance().logLevel) {
            std::string levelStr;
            switch (level) {
                case LogLevel::DEBUG:
                    levelStr = "DEBUG";
                    break;
                case LogLevel::INFO:
                    levelStr = "INFO";
                    break;
                case LogLevel::WARNING:
                    levelStr = "WARNING";
                    break;
                case LogLevel::ERROR:
                    levelStr = "ERROR";
                    break;
            }
            seqan3::debug_stream << "[" << levelStr << "] " << getTime() << " ";
            (seqan3::debug_stream << ... << std::forward<Args>(args)) << "\n";

            if (level == LogLevel::ERROR) {
                exit(EXIT_FAILURE);
            }
        }
    }

   private:
    Logger() : logLevel(LogLevel::INFO) {}
    Logger(const Logger &) = delete;
    Logger &operator=(const Logger &) = delete;

    LogLevel logLevel;
    std::mutex logMutex;

    static std::string getTime() {
        const auto now = std::chrono::system_clock::now();
        const std::time_t current_time = std::chrono::system_clock::to_time_t(now);

        std::ostringstream time_stream;
        time_stream << std::put_time(std::localtime(&current_time), "[%Y-%m-%d %H:%M:%S]") << " ";

        return time_stream.str();
    };
};