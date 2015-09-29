#pragma once

#include <string>
#include <iostream>

class Logger {
public:
  static Logger& get_instance() {
    static Logger instance;

    return instance;
  }

  void debug(std::string _message);
  void info(std::string _message);
  void warn(std::string _message);
  void error(std::string _message);
  void fatal(std::string _message);

private:
  Logger() {}

  Logger(Logger const&) = delete;
  void operator=(Logger const&) = delete;
};
