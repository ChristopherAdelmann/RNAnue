#include "Helper.hpp"

#include <iostream>
#include <string>

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
    static Logger &getInstance()
    {
        static Logger instance; // Singleton instance
        return instance;
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