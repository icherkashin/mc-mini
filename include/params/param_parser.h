#pragma once

#include <vector>
#include <string>

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>

#include "params.h"
#include "params/param_node.h"

class ParamParser {
public:
  ParamParser();

  static Params& parse(std::string filename);

private:
  static void handle_line(ParamNode *(&cursor), std::string line);

  static void enter_action(ParamNode *(&cursor), std::string rest);
  static void leave_action(ParamNode *(&cursor), std::string rest);
  static void set_action(ParamNode *(&cursor), std::string rest);
};
