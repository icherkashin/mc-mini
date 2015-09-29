#include <fstream>
#include <iostream>

#include <vector>
#include <string>

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>

#include "logger.h"
#include "params/param_parser.h"

static Logger &Logger = Logger::get_instance();

/**
 * Parse a file for parameters.
 */
Params& ParamParser::parse(std::string _filename) {
  // Open the file for parsing
  std::ifstream paramfile(_filename);

  logger.debug("Parsing parameter file '" + _filename + "'");

  // Use a cursor to represent the current position in the parameter tree
  ParamNode *cursor = new ParamNode();

  // Read the file line by line
  std::string line;
  while (getline(paramfile, line)) handle_line(cursor, line);

  // Close the file once finished
  paramfile.close();

  logger.debug("Finished parsing parameters.");

  return (*(new Params(cursor)));
}

/**
 * Parse a single line of a parameter file
 */
void ParamParser::handle_line(ParamNode *(&_cursor), std::string _line) {
  // Find the first occurence of a space in the string
  boost::algorithm::trim(_line);
  const size_t divider_idx = _line.find_first_of(" \t\n");

  std::string action;
  std::string rest;

  // Divide out the 'action' and the 'value', based ont he position of
  // the first divider.
  if (divider_idx != std::string::npos) {
    action = _line.substr(0, divider_idx);
    rest = _line.substr(divider_idx + 1);
  } else {
    action = _line;
  }

  // strip any comments from the remainder of the string and trim whitespace
  const size_t comment_idx = rest.find_first_of("#");
  if (comment_idx != std::string::npos) {
    rest = rest.substr(0, comment_idx);
  }
  boost::algorithm::trim(rest);

  // We can ignore the line if it was either empty or a comment
  if (action == "" || action[0] == '#') return;

  // Otherwise, handle the action
  if (action == "enter") {
    enter_action(_cursor, rest);
  } else if (action == "leave") {
    leave_action(_cursor, rest);
  } else if (action == "set") {
    set_action(_cursor, rest);
  } else {
    logger.error("Couldn't parse param file");
    // We got something that we can't parse...
    throw("error");
  }
}

void ParamParser::enter_action(
        ParamNode *(&_cursor),
        std::string _rest) {

  // parse out the section title
  std::string section_key;
  std::string rest;

  const size_t divider_idx = _rest.find_first_of(" \t\n");

  if (divider_idx != std::string::npos) {
    section_key = _rest.substr(0, divider_idx);
    rest = _rest.substr(divider_idx + 1);
  } else {
    section_key = _rest;
  }

  // The remainder of the line should be empty or a comment
  boost::algorithm::trim(rest);
  if (rest != "") throw("error");

  ParamNode *section = _cursor->get_child(section_key);
  // Make sure that the required section exists
  if (section == NULL) {
    logger.debug("Creating new parameter section '" + section_key + "'");

    // If the section doesn't exist, create it.
    ParamNode *new_section = ParamNode(_cursor);
    _cursor->set_child(section_key, new_section);
    section = new_section;
  }

  logger.debug("Entering parameter section '" + section_key + "'");

  // Move the cursor to the entered section.
  _cursor = section;
}

void ParamParser::leave_action(
        ParamNode *(&_cursor),
        std::string _rest) {
  // The remainder of the string should be empty
  if (_rest != "") throw("error");

  // We shouldn't be able to leave the root node
  if (_cursor = NULL) throw("error");

  logger.debug("Leaving parameter section");

  // Leave the current section
  _cursor = _cursor->get_parent();
}

void ParamParser::set_action(
        ParamNode *(&_cursor),
        std::string _rest) {

  std::string data_key;
  std::string data;

  // Find the position of the '=' sign
  const size_t equals_idx = _rest.find_first_of("=");
  if (equals_idx != std::string::npos) {
    data_key = _rest.substr(0, equals_idx);
    _rest = _rest.substr(equals_idx + 1);
  } else {
    data_key = _rest;
  }

  boost::algorithm::trim(data_key);
  boost::algorithm::trim(data);

  // There must be a data key and data entry
  if (data_key == "" || data_entry == "") throw("error");

  logger.debug("Setting the value of parameter '" + data_key + "' to '" + data + "'");

  _cursor->set_data(data_key, data);
}
