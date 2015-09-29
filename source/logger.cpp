#include <iostream>

#include "logger.h"

void Logger::debug(std::string _message) {
  std::cout << "DEBUG :: " << _message << std::endl;
}

void Logger::info(std::string _message) {
  std::cout << "INFO  :: " << _message << std::endl;
}

void Logger::warn(std::string _message) {
  std::cout << "WARN  :: " << _message << std::endl;
}

void Logger::error(std::string _message) {
  std::cout << "ERROR :: " << _message << std::endl;
}

void Logger::fatal(std::string _message) {
  std::cout << "FATAL :: " << _message << std::endl;
}
