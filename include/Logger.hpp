#ifndef LOGGER_H
#define LOGGER_H

#include "Helper.hpp"

#include <iostream>
#include <string>
#include <map>
#include <optional>

#include <seqan3/core/debug_stream.hpp>

enum class LogLevel
{
    INFO,
    WARNING,
    ERROR
};

class Logger
{
public:
    static std::optional<LogLevel> stringToLogLevel(const std::string &level)
    {
        static const std::map<std::string, LogLevel> stringToLogLevelMap{
            {"info", LogLevel::INFO},
            {"warning", LogLevel::WARNING},
            {"error", LogLevel::ERROR}};

        auto it = stringToLogLevelMap.find(level);
        if (it == stringToLogLevelMap.end())
        {
            return std::nullopt;
        }
        return it->second;
    }

    static Logger &getInstance()
    {
        static Logger instance; // Singleton instance
        return instance;
    }

    static void setLogLevel(const std::string &logLevelString)
    {
        std::optional<LogLevel> logLevel = stringToLogLevel(logLevelString);
        if (logLevel.has_value())
        {
            getInstance().i_setLogLevel(logLevel.value());
        }
        else
        {
            Logger::log(LogLevel::ERROR, "Invalid log level: " + logLevelString);
            exit(1);
        }
    }

    // Convenience functions for setting log level and logging without calling getInstance
    static void setLogLevel(LogLevel level)
    {
        getInstance().i_setLogLevel(level);
    }

    template <typename T>
    static void log(LogLevel level, const T &message)
    {
        getInstance().i_log(level, message);
    }

private:
    Logger() : logLevel(LogLevel::INFO) {}      // Private constructor to prevent instantiation
    Logger(const Logger &) = delete;            // Delete copy constructor
    Logger &operator=(const Logger &) = delete; // Delete assignment operator

    LogLevel logLevel;

    void i_setLogLevel(LogLevel level)
    {
        logLevel = level;
    }

    template <typename T>
    void i_log(LogLevel level, const T &message)
    {
        if (level >= logLevel)
        {
            std::string levelStr;
            switch (level)
            {
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
            seqan3::debug_stream << "[" << levelStr << "] " << helper::getTime() << " " << message << std::endl;
        }
    }
};

#endif // LOGGER_H