#pragma once

#include <exception>

using namespace std;

class param_not_found : public exception
{
  virtual const char * what() const throw() {
    return "Parameter not found";
  }
} param_not_found;

class child_not_found : public exception
{
  virtual const char * what() const throw() {
    return "Child not found";
  }
} child_not_found;
