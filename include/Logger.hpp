#pragma once

// Standard
#include <chrono>
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
            {"info", LogLevel::INFO}, {"warning", LogLevel::WARNING}, {"error", LogLevel::ERROR}};

        auto it = stringToLogLevelMap.find(logLevelString);
        if (it == stringToLogLevelMap.end()) {
            Logger::log(LogLevel::ERROR, "Invalid log level: " + logLevelString);
            exit(1);
        }
        getInstance().logLevel = it->second;
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
        }
    }

   private:
    Logger() : logLevel(LogLevel::INFO) {}       // Private constructor to prevent instantiation
    Logger(const Logger &) = delete;             // Delete copy constructor
    Logger &operator=(const Logger &) = delete;  // Delete assignment operator

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